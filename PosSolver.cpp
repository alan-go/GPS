#include "PosSolver.h"
#include "GNSS.h"

PosSolver::PosSolver(SVs svs, char *raw, NtripRTK *rtk, GNSS *gnss): svs(svs),raw(raw),rtk(rtk),gnss(gnss){
    xyz = gnss->xyz;
    tu = gnss->tu;
    tuBeiDou = gnss->tuBeiDou;
    tuGps = gnss->tuGps;
}

int PosSolver::CalcuPosition() {
    PrepareData(svs, raw);

    int countBeiDou = 0, countGPS = 0;
    for(int i=0;i<numMeas;i++)
    {
        SV *svTemp = visibleSvs[i];
        printf("svsfor visiable:%d,%02d:%d,%d,%d.health:%d\n",
               svTemp->type,svTemp->svId,svTemp->page1OK,svTemp->page2OK,svTemp->page3OK,svTemp->SatH1);
        if(svTemp->JudgeUsable(gnss->useBeiDou,gnss->useGPS)){
            SvsForCalcu.push_back(svTemp);
            if(SV::BeiDou == svTemp->type)countBeiDou++;
            if(SV::GPS == svTemp->type)countGPS++;
            printf("svsfor calcu\n");
        }
    }
    if(SvsForCalcu.size()<4){
        printf("calcu:Not enough Svs.\n");
        return -1;
    } else if(4==SvsForCalcu.size()&&countBeiDou*countBeiDou){
        printf("4 SVs with GPS and BeiDou,\nUnable to solve.\n");
    } else if(0==countBeiDou*countGPS){
        solvePosition();
    } else{
        solvePositionBeiDouGPS();
    }
}

int PosSolver::PrepareData(SVs svs, char *raw) {
    char* playload = raw + 6;

    char* temp = playload;
    rcvtow = *(double*)temp;
    temp = playload + 11;
    numMeas = *(u_int8_t*)temp;
    if(0==numMeas)return 1;

    for(u_int8_t n = 0;n<numMeas;n++){
        SV* svTemp;
        int n32 = n*32;
        uint8_t gnssId = *(uint8_t *)(playload+36+n32);
        uint8_t svId = *(uint8_t *)(playload+37+n32);
        //watchable SVs
        switch (gnssId){
            case 0:
                svTemp = &(svs.svGpss[svId-1]);
                break;
            case 3:
                svTemp = &(svs.svBeiDous[svId-1]);
                break;
        }
        visibleSvs.push_back(svTemp);
        svTemp->prMes = *(double*)(playload+16+n32);
        svTemp->cpMes = *(double*)(playload+24+n32);
        svTemp->doMes = *(float *)(playload+32+n32);

    }
    printf("start solving rawdata numMesa=%d\n",numMeas);
    return 0;
}

int PosSolver::solvePosition() {
    int N = SvsForCalcu.size();
    Vector4d dtxyzt = Vector4d::Ones();
    MatrixXd pc(N,1);
    MatrixXd b(N,1);
    cout<<"-----------------start clacu svPosition\n-------------"<<endl;

    for(int i = 0;i< N;i++){
        SV *sv= SvsForCalcu[i];
        sv->CalcuECEF(rcvtow);
        //todo: calcu I,T,dtu
        double pci = sv->prMes + Light_speed * sv->tsDelta - sv->I - sv->T;
        pc(i) = pci;
        sv->PrintInfo(1);
        cout<<"norm"<<sv->position.norm()<<endl;

    }
    cout<<"+++++++++++++ poss\n-------------"<<endl;

    int numCalcu = 0;
    while (dtxyzt.norm()>0.1){
        numCalcu++;
        MatrixXd G(N,4);
        for(int i = 0;i< N;i++){
            SV *sv= SvsForCalcu[i];
            double r = (sv->position-xyz).norm();
            //自转修正,这个我需要考虑一下是不是在这里处理
            double omegat = r/Light_speed*Omega_e;
            MatrixXd earthRotate(3,3);
            earthRotate<<cos(omegat),sin(omegat),0,-sin(omegat),cos(omegat),0,0,0,1;
            Vector3d svPositionEarthRotate = earthRotate*sv->position;
            cout<<"positon = \n"<<sv->position<<endl;
            cout<<"earthRotate = \n"<<earthRotate<<endl;
            cout<<"positionRotate = \n"<<svPositionEarthRotate<<endl;
            r = (svPositionEarthRotate-xyz).norm();
//            printf("r=%.10f\n",r);
//            r = (sv->position-rxyz).norm()+(sv->position(0)*rxyz(1)-sv->position(1)*rxyz(0))*Omega_e/Light_speed;
//            printf("r=%.10f\n",r);
//            printf("r=%.10f\n",(sv->position-rxyz).norm());

            b(i) = pc(i)-r-tu;  //这里的r需要考虑一下
            G(i,0) = (xyz(0)-svPositionEarthRotate(0))/r;  //x
            G(i,1) = (xyz(1)-svPositionEarthRotate(1))/r;  //y
            G(i,2) = (xyz(2)-svPositionEarthRotate(2))/r;  //z
            G(i,3) = 1;
        }
        MatrixXd GT = G.transpose();
        dtxyzt = ((GT*G).inverse())*GT*b;

        cout<<"b = "<<b<<endl;
        cout<<"dt = "<<dtxyzt<<endl;
        xyz += dtxyzt.head(3);
        tu+=dtxyzt(3);
        cout<<"++++++++++++xyzt\n"<<xyz<<endl;
        cout<<"xyzNorm\n"<<xyz.norm()<<endl;
        XYZ2LLA(xyz,LLA);
        cout<<"++++++++++++LLA\n"<<LLA<<endl;
        if(numCalcu>30)return false;

    }

    XYZ2LLA(xyz.head(3),LLA);
    cout<<"++++++++++++rxyzt\n"<<xyz<<endl;
    cout<<"++++++++++++LLA\n"<<LLA<<endl;
    return true;

}

int PosSolver::solvePositionBeiDouGPS(){
    int N = SvsForCalcu.size();
    VectorXd dtxyzBG(5);
    dtxyzBG<<1,1,1,1,1;
    MatrixXd pc(N,1);
    MatrixXd b(N,1);
    cout<<"-----------------start clacu svPosition\n-------------"<<endl;

    for(int i = 0;i< N;i++){
        SV *sv= SvsForCalcu[i];
        sv->CalcuECEF(rcvtow);
        //todo: calcu I,T,dtu
        double pci = sv->prMes + Light_speed * sv->tsDelta - sv->I - sv->T;
        pc(i) = pci;
        sv->PrintInfo(1);
        cout<<"norm"<<sv->position.norm()<<endl;
    }
    cout<<"+++++++++++++ poss\n-------------"<<endl;

    int numCalcu = 0;
    while (dtxyzBG.norm()>0.1){
        numCalcu++;
        MatrixXd G(N,5);
        for(int i = 0;i< N;i++){
            SV *sv= SvsForCalcu[i];
            cout<<"rxyz=\n"<<xyz<<endl;
            double r = (sv->position-xyz).norm();
            //自转修正,这个我需要考虑一下是不是在这里处理
            double omegat = r/Light_speed*Omega_e;
            Matrix<double,3,3>earthRotate;
            earthRotate<<cos(omegat),sin(omegat),0,-sin(omegat),cos(omegat),0,0,0,1;
            Vector3d svPositionEarthRotate = earthRotate*sv->position;
            cout<<"positon = \n"<<sv->position<<endl;
            cout<<"earthRotate = \n"<<earthRotate<<endl;
            cout<<"positionRotate = \n"<<svPositionEarthRotate<<endl;
            r = (svPositionEarthRotate-xyz).norm();
//            printf("r=%.10f\n",r);
//            r = (sv->position-rxyz).norm()+(sv->position(0)*rxyz(1)-sv->position(1)*rxyz(0))*Omega_e/Light_speed;
//            printf("r=%.10f\n",r);
//            printf("r=%.10f\n",(sv->position-rxyz).norm());

            b(i) = pc(i)-r-tuBeiDou-tuGps;  //这里的r需要考虑一下
            G(i,0) = (xyz(0)-svPositionEarthRotate(0))/r;  //x
            G(i,1) = (xyz(1)-svPositionEarthRotate(1))/r;  //y
            G(i,2) = (xyz(2)-svPositionEarthRotate(2))/r;  //z
            G(i,3) = (SV::BeiDou==sv->type)?1:0;
            G(i,4) = (SV::GPS==sv->type)?1:0;
        }
        MatrixXd GT = G.transpose();
        dtxyzBG = ((GT*G).inverse())*GT*b;

        cout<<"G=\n"<<G<<endl;
        cout<<"GTG=\n"<<GT*G<<endl;
        cout<<"GTGi=\n"<<(GT*G).inverse()<<endl;
        cout<<"b = "<<b<<endl;
        cout<<"dt = "<<dtxyzBG<<endl;
        xyz += dtxyzBG.head(3);
        tuBeiDou+=dtxyzBG(3);
        tuGps+=dtxyzBG(4);
        cout<<"++++++++++++xyzt\n"<<xyz<<endl;
        cout<<"xyzNorm\n"<<xyz.norm()<<endl;
        XYZ2LLA(xyz,LLA);
        cout<<"++++++++++++LLA\n"<<LLA<<endl;
        if(numCalcu>30)return false;

    }

    XYZ2LLA(xyz,LLA);
    cout<<"++++++++++++rxyzt\n"<<xyz<<endl;
    cout<<"++++++++++++LLA\n"<<LLA<<endl;
    return true;
}

int PosSolver::XYZ2LLA(Vector3d XYZ, Vector3d &LLA) {
    double x,y,z,r,sinA,cosA,sinB,cosB,N,h,lati;
    int i;
    double v_xyz[3],v_enu[3],a_xyz[3],a_enu[3];
    double lati_f,long_f;   // N_geoi表中的经纬度，以角度记，均加到正值处以10
    double lati_w,long_w;   //权重
    int lati_index,long_index,long_index_n;  //用于N_geoid中的下标
    //input
    x = XYZ(0);
    y = XYZ(1);
    z = XYZ(2);
    v_xyz[0] = v_xyz[0];
    v_xyz[1] = v_xyz[1];
    v_xyz[2] = v_xyz[2];
    a_xyz[0] = a_xyz[0];
    a_xyz[1] = a_xyz[1];
    a_xyz[2] = a_xyz[2];
    //起始迭代点
    r = sqrt(x*x+y*y);
    h = 0;
    lati = 0;
    for (i=0;i<6;i++){
        sinA = sin(lati);
        cosA = cos(lati);
        N = Earth_a / sqrt(1-sinA*sinA*Earth_ee);
        h = r/cosA - N;
        lati = atan(z/(r*(1-Earth_ee*N/(N+h))));
    }
    //output Longtitude Latitude High  SS[8] 转换矩阵
    LLA(0) = atan2(y,x)*180/GPS_PI;
    LLA(1) = lati*180/GPS_PI;
    LLA(2) = h;

    //给转换矩阵赋值
    sinB = y/r;
    cosB = x/r;
    double SS[8];
    SS[0] = -sinB;
    SS[1] = cosB;
    SS[2] = 0;
    SS[3] = -sinA*cosB;
    SS[4] = -sinA*sinB;
    SS[5] = cosA;
    SS[6] = cosA*cosB;
    SS[7] = cosA*sinB;
    SS[8] = sinA;

    //XYZ向速度，加速度ENU 之后有时间自己实现吧，下面就是一个矩阵乘法然后赋值
    //SS*v_xyz = v_enu   SS*a_xyz = a_enu
    return 0;
}