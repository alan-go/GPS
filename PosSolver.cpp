#include "PosSolver.h"
#include "GNSS.h"

PosSolver::PosSolver(SVs svs, NtripRTK *rtk, GNSS *gnss): svs(svs),rtk(rtk),gnss(gnss),numBDSUsed(0),numGPSUsed(0){
    if(gnss->isPositioned){
        xyz = gnss->records[0].xyz;
        LLA = gnss->records[0].lla;
    }
    tu = gnss->tu;
    tuBeiDou = gnss->tuBeiDou;
    tuGps = gnss->tuGps;
}

PosSolver::~PosSolver(){}

int PosSolver::PrepareSVsData(vector<SV *> &svsForCalcu) {
    ReadVisibalSvsRaw(svs, visibleSvs, raw);
    numGPSUsed = 0;numBDSUsed = 0;
    SelectSvsFromVisible(visibleSvs,svsForCalcu);

    for(int i = 0;i< svsForCalcu.size();i++){
        int timeForECEF = rcvtow;
        SV *sv= svsForCalcu[i];
        if(SV::BeiDou == sv->type)timeForECEF-=14;
        sv->CalcuTime(timeForECEF);
        sv->CalcuECEF(sv->tsReal);
        if(gnss->isPositioned)sv->CorrectIT(xyz,LLA,timeForECEF);

        XYZ2LLA(sv->position,sv->sLLA);
//        sv->PrintInfo(2);
//        sv->PrintInfo(1);
//        cout<<"norm"<<sv->position.norm()<<endl;
    }
}

int PosSolver::PositionRtk() {
    int sigInd = 1;
    //插值
    int result = 0;
    vector<SV*> svs0, svs1;
    PrepareSVsData(svs0);
    printf("rcvtow = %lf\n",rcvtow);
    //暂时单模计算
    if(numGPSUsed*numBDSUsed){
        printf("-------------------count 0\n");
        return -1;
    }
    SV* svCectre;
    int svCentreInd;
    double elevMax = 0;
//    gnss->rtkManager.mtxData.lock();
    printf(" positon rtk lock 0\n");
    printf("---------calcu N0 = %d\n",svs0.size());

    for(int i = 0,k=-1; i<svs0.size();i++){
        SV* sv = svs0[i];
        int time = rcvtow;
        if(sv->type==SV::BeiDou)time-=14;
        double pr = sv->InterpRtkData(time,sigInd);
        if(0!=pr){
            k++;
            svs1.push_back(sv);
            if(sv->elevationAngle>elevMax){
                elevMax = sv->elevationAngle;
                svCectre = sv;
                svCentreInd = k;
            }
        }
    }
    printf("---------calcu N = %d\n",svs1.size());
    svs1.erase(svs1.begin()+svCentreInd);

    int N = svs1.size();
    printf("---------calcu N = %d\n",N);
    if(N<3){
        gnss->rtkManager.mtxData.unlock();
        return -1;
    }
    MatrixXd pur(N,1);
    MatrixXd Gur(N,3);

    double pur0 = svCectre->prMes - svCectre->prInterp[sigInd];
    Vector3d referBase = gnss->rtkManager.ECEF_XYZ;
    Vector3d Ir0 = svCectre->position - referBase;
    Ir0 = Ir0/Ir0.norm();

    for(int i=0;i<svs1.size();i++){
        SV* sv = svs1[i];
        pur(i) = (sv->prMes-sv->prInterp[sigInd])-pur0;
        Vector3d Iri = sv->position - referBase;
        Iri = Iri/Iri.norm();
        Vector3d Ir_0i = Ir0 - Iri;
        Gur(i,0) = Ir_0i(0);
        Gur(i,1) = Ir_0i(1);
        Gur(i,2) = Ir_0i(2);
    }
    printf(" positon rtk unlock 0\n");
//    gnss->rtkManager.mtxData.unlock();
    printf(" positon rtk unlock 1\n");

    MatrixXd GurT = Gur.transpose();
    Vector3d bur = (GurT*Gur).inverse()*GurT*pur;
    xyz = referBase+bur;


    XYZ2LLA(xyz,LLA);
    gnss->isPositioned = true;
    gnss->AddPosRecord(GNSS::PosRcd(rcvtow,xyz,LLA));
    cout<<"++++++++++++xyz\n"<<xyz<<endl;
    printf("++++++++++++LLA === %lf,%lf,%lf\n",LLA(1)*180/GPS_PI,LLA(0)*180/GPS_PI,LLA(2));
}

int PosSolver::PositionSingle() {
    int result = 0;
    vector<SV*> svsForCalcu;
    PrepareSVsData(svsForCalcu);
    printf("rcvtow = %lf\n",rcvtow);
    int N = svsForCalcu.size();
    if(N<4){
        printf("\n\n----------calcu:%d,Not enough Svs.\n",N);
        return -1;
    } else if(4==N && numGPSUsed*numBDSUsed){
        printf("4 SVs with GPS and BeiDou:%d, %d,   Unable to solve.\n",numGPSUsed,numBDSUsed);
        return -1;
    } else if(0==numBDSUsed*numGPSUsed){
        result = SolvePosition(svsForCalcu);
    } else{
        result = SolvePositionBeiDouGPS(svsForCalcu);
    }

//    if(result)
    delete(this);
}


int PosSolver::SolvePosition(vector<SV*>svsForCalcu) {
    int N = svsForCalcu.size();
    printf("\nSolve Position N = %d\n",N);
    Vector4d dtxyzt = Vector4d::Ones();
    MatrixXd pc(N,1);
    MatrixXd b(N,1);

    for(int i = 0;i< N;i++){
        SV *sv= svsForCalcu[i];
        pc(i) = sv->prMes + Light_speed * sv->tsDelta - sv->I - sv->T;
    }

    int numCalcu = 0;
    MatrixXd G(N,4);
    Matrix<double,4,4> H;
    while (dtxyzt.norm()>0.1){
        numCalcu++;
        for(int i = 0;i< N;i++){
            SV *sv= svsForCalcu[i];
            double r = (sv->position-xyz).norm();
            //todo 自转修正
            double omegat = -r/Light_speed*Omega_e;
//            printf("r=%lf,omegat = %lf\n",r,omegat);
//            double omegat = 0;
            MatrixXd earthRotate(3,3);
            earthRotate<<cos(omegat),sin(omegat),0,-sin(omegat),cos(omegat),0,0,0,1;
            Vector3d svPositionEarthRotate = earthRotate*sv->position;
//            cout<<"positon = \n"<<sv->position<<endl;
//            cout<<"earthRotate = \n"<<earthRotate<<endl;
//            cout<<"positionRotate = \n"<<svPositionEarthRotate<<endl;
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
        H = (GT*G).inverse();
        dtxyzt = H*GT*b;

//        cout<<"b = "<<b<<endl;
        cout<<"dt = "<<dtxyzt<<endl;
        xyz += dtxyzt.head(3);
        tu+=dtxyzt(3);
//        cout<<"++++++++++++xyzt\n"<<xyz<<endl;
//        cout<<"xyzNorm\n"<<xyz.norm()<<endl;
//        XYZ2LLA(xyz,LLA);
//        cout<<"++++++++++++LLA\n"<<LLA<<endl;
        printf("num = %d, norm = %lf\n",numCalcu,dtxyzt.norm());
        if(numCalcu>30)return false;
    }
    double PDOP = sqrt(H(0,0)+H(1,1)+H(2,2));
    XYZ2LLA(xyz,LLA);
    gnss->isPositioned = true;
    gnss->AddPosRecord(GNSS::PosRcd(rcvtow,xyz,LLA));
    gnss->tu = tu;
    cout<<"++++++++++++xyz\n"<<xyz<<endl;
    printf("++++++++++++LLA === %lf,%lf,%lf\n",LLA(1)*180/GPS_PI,LLA(0)*180/GPS_PI,LLA(2));
//    cout<<"++++++++++++LLA\n"<<LLA*180/GPS_PI<<endl;

    printf("rcvtow = %lf\n",rcvtow);
    for(int i=0;i<N;i++)
    {
        double prres = b(i);
        SV *svTemp = svsForCalcu[i];
//        svTemp->CalcuelEvationAzimuth(xyz,LLA);
        printf("svs visiable:%d,%02d:pr:%lf\tprres:%.4f,\t%f\n",
               svTemp->type,svTemp->svId,svTemp->prMes,prres,svTemp->elevationAngle);
    }
    
    printf("PDOP = %lf\n",PDOP);
    return 1;
}

int PosSolver::SolvePositionBeiDouGPS(vector<SV*>svsForCalcu){
    int N = svsForCalcu.size();
    printf("\nSolve Position GPS-BeiDou, N = %d\n",N);
    VectorXd dtxyzBG(5);
    Matrix<double, 5,5> H;
    dtxyzBG<<1,1,1,1,1;
    MatrixXd pc(N,1);
    MatrixXd b(N,1);

    for(int i = 0;i< N;i++){
        SV *sv= svsForCalcu[i];
        pc(i) = sv->prMes + Light_speed * sv->tsDelta - sv->I - sv->T;
    }

    int numCalcu = 0;
    while (dtxyzBG.norm()>0.1){
        numCalcu++;
        MatrixXd G(N,5);
        for(int i = 0;i< N;i++){
            SV *sv= svsForCalcu[i];
//            cout<<"rxyz=\n"<<xyz<<endl;
            double r = (sv->position-xyz).norm();
            //自转修正,这个我需要考虑一下是不是在这里处理
            double omegat = r/Light_speed*Omega_e;
            Matrix<double,3,3>earthRotate;
            earthRotate<<cos(omegat),sin(omegat),0,-sin(omegat),cos(omegat),0,0,0,1;
            Vector3d svPositionEarthRotate = earthRotate*sv->position;
//            cout<<"positon = \n"<<sv->position<<endl;
//            cout<<"earthRotate = \n"<<earthRotate<<endl;
//            cout<<"positionRotate = \n"<<svPositionEarthRotate<<endl;
            r = (svPositionEarthRotate-xyz).norm();
//            printf("r=%.10f\n",r);
//            r = (sv->position-rxyz).norm()+(sv->position(0)*rxyz(1)-sv->position(1)*rxyz(0))*Omega_e/Light_speed;
//            printf("r=%.10f\n",r);
//            printf("r=%.10f\n",(sv->position-rxyz).norm());

            G(i,0) = (xyz(0)-svPositionEarthRotate(0))/r;  //x
            G(i,1) = (xyz(1)-svPositionEarthRotate(1))/r;  //y
            G(i,2) = (xyz(2)-svPositionEarthRotate(2))/r;  //z
            G(i,3) = (SV::BeiDou==sv->type)?1:0;
            G(i,4) = (SV::GPS==sv->type)?1:0;
//            b(i) = pc(i)-r-tuBeiDou-tuGps;  //这里的r需要考虑一下
            b(i) = pc(i)-r-tuBeiDou*G(i,3)-tuGps*G(i,4);  //这里的r需要考虑一下

        }
        MatrixXd GT = G.transpose();
        H = (GT*G).inverse();
        dtxyzBG = H*GT*b;

//        cout<<"G=\n"<<G<<endl;
//        cout<<"GTG=\n"<<GT*G<<endl;
//        cout<<"GTGi=\n"<<(GT*G).inverse()<<endl;
//        cout<<"b = "<<b<<endl;
        cout<<"dt = "<<dtxyzBG<<endl;
        xyz += dtxyzBG.head(3);
        tuBeiDou+=dtxyzBG(3);
        tuGps+=dtxyzBG(4);
//        cout<<"++++++++++++xyzt\n"<<xyz<<endl;
//        cout<<"xyzNorm\n"<<xyz.norm()<<endl;
//        XYZ2LLA(xyz,LLA);
//        cout<<"++++++++++++LLA\n"<<LLA<<endl;
        printf("num = %d, norm = %lf\n",numCalcu,dtxyzBG.norm());
        if(numCalcu>30)return false;
    }
    XYZ2LLA(xyz,LLA);
    gnss->isPositioned = true;
    gnss->AddPosRecord(GNSS::PosRcd(rcvtow,xyz,LLA));
    gnss->tuGps = tuGps;
    gnss->tuBeiDou = tuBeiDou;
    cout<<"++++++++++++rxyzt\n"<<xyz<<endl;
    printf("++++++++++++LLA === %lf,%lf,%lf\n",LLA(1)*180/GPS_PI,LLA(0)*180/GPS_PI,LLA(2));
    double PDOP = sqrt(H(0,0)+H(1,1)+H(2,2));
    printf("PDOP = %lf\n",PDOP);

    printf("rcvtow = %lf\n",rcvtow);
    for(int i=0;i<N;i++)
    {
        double prres = b(i);
        SV *svTemp = svsForCalcu[i];
//        svTemp->CalcuelEvationAzimuth(xyz,LLA);
        printf("svs visiable:%d,%02d:pr:%lf\tprres:%.4f,\t%f\n",
               svTemp->type,svTemp->svId,svTemp->prMes,prres,svTemp->elevationAngle);
    }

    return 1;
}



int PosSolver::SolvePositionCalman() {

}



int PosSolver::SelectSvsFromVisible(vector<SV*> &all, vector<SV*> &select) {
    //first judge ephemeric
    for(int i=0;i<numMeas;i++)
    {
        SV *svTemp = all[i];

        printf("svs visiable:%d,%02d:%d,%d,%d.health:%d,pr=%lf\n",
               svTemp->type,svTemp->svId,svTemp->bstEphemOK[0],svTemp->bstEphemOK[1],
               svTemp->bstEphemOK[2],svTemp->SatH1,svTemp->prMes);
        bool useForCalcu = svTemp->JudgeUsable(gnss->useBeiDou,gnss->useGPS);
        useForCalcu = useForCalcu && svTemp->MeasureGood();
        if(gnss->isPositioned)useForCalcu = useForCalcu && svTemp->ElevGood();
        if(svTemp->AODC>24){
            printf("\n\n\nAAAAAODC = %d\n\n\n",svTemp->AODC);
            useForCalcu = 0;
        }
        if(useForCalcu){
            select.push_back(svTemp);
            if(SV::BeiDou == svTemp->type)numBDSUsed++;
            if(SV::GPS == svTemp->type)numGPSUsed++;
            printf("svsfor calcu\n");
        }
    }
}


int PosSolver::ReadVisibalSvsRaw(SVs svs,vector<SV*> &svVisable, char *raw) {
    char* playload = raw + 6;

    char* temp = playload;
    rcvtow = *(double*)temp;
    temp = playload + 11;
    numMeas = *(u_int8_t*)temp;
    printf("prepare rawdata , numMesa=%d\n",numMeas);
    if(0==numMeas)return -1;

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
        svVisable.push_back(svTemp);
        svTemp->prMes = *(double*)(playload+16+n32);
        svTemp->cpMes = *(double*)(playload+24+n32);
        svTemp->doMes = *(float *)(playload+32+n32);

    }
    return 0;
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
    for (i=0;i<5;i++){
        sinA = sin(lati);
        cosA = cos(lati);
        N = Earth_a / sqrt(1-sinA*sinA*Earth_ee);
        h = r/cosA - N;
        lati = atan(z/(r*(1-Earth_ee*N/(N+h))));
//        printf("lati = %lf\n",lati);
    }
    //output Longtitude Latitude High  SS[8] 转换矩阵
    LLA(0) = atan2(y,x);
    LLA(1) = lati;
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