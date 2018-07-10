#include "UbloxSolver.h"

using namespace std;

SvInfo::SvInfo(){
    page1OK = page2OK = page3OK = false;
    I = T = 0;
}
SvInfo::~SvInfo(){}

bool SvInfo::CalcuTime(double rcvtow) {
    ts = rcvtow - prMes/Light_speed;
    double t = ts-toc;
    //此时的钟差是没有考虑相对论效应和 TGD的
    for(int i = 0;i<5;i++){
        tsDelta = a0+a1*t+a2*t*t;
        t -=tsDelta;
    }
    tsReal = ts - tsDelta;
    return true;
}

void SvInfo::PrintInfo(int printType) {
    if(1==printType){
        cout<<"+++++++SvPosition:"<<type<<","<<svId<<endl;
        cout<<position<<endl;
        cout<<"dt"<<tsDelta<<endl;
    }
}


 bool SvInfo::CalcuECEF(double rcvtow) {
    CalcuTime(rcvtow);
    double A = orbit.sq_a*orbit.sq_a;
    double n0 = sqrt(M_miu/(A*A*A));
    double tk = tsReal - orbit.toe;
    if(tk > 302400)tk-=604800;
    if(tk < -302400)tk+=604800;
    double n = n0 + orbit.dtn;
    double Mk = orbit.M0 + n*tk;
    double Ek = Mk,EkOld = Ek-1;
    while(abs(Ek-EkOld)>1e-8){
        EkOld = Ek;
        Ek = EkOld-(EkOld-orbit.e*sin(EkOld)-Mk)/(1-orbit.e*cos(EkOld));
    }
    cout<<"EK:"<<Ek<<endl;
    double vk = atan((sqrt(1-orbit.e*orbit.e)*sin(Ek))/(cos(Ek)-orbit.e));
    double phyk = vk + orbit.omega;
    double sin2phy = sin(2*phyk),cos2phy = cos(2*phyk);
    double dtuk = orbit.Cus*sin2phy + orbit.Cuc*cos2phy;
    double dtrk = orbit.Crs*sin2phy + orbit.Crc*cos2phy;
    double dtik = orbit.Cis*sin2phy + orbit.Cic*cos2phy;
    double uk = phyk + dtuk;
    double rk = A*(1-orbit.e*cos(Ek)) + dtrk;
    double ik = orbit.i0 + orbit.IDOT*tk + dtik;
    double xk = rk * cos(uk);
    double yk = rk * sin(uk);
    double Omegak = orbit.Omega0 + (orbit.OmegaDot - isBeiDouGEO?0:Omega_e) * tk - Omega_e * orbit.toe;
    MatrixXd transfer(3,2);
    transfer<<cos(Omegak),-cos(ik)*sin(Omegak),sin(Omegak),cos(ik)*cos(Omegak),0,sin(ik);
    if(isBeiDouGEO){
        Vector3d xyzGK = transfer*Vector2d(xk,yk);
        double phyX = -5/180*M_PI;
        double phyZ = Omega_e * tk;
        Matrix3d Rz,Rx;
        Rx<<1,0,0,0,cos(phyX),sin(phyX),0,-sin(phyX),cos(phyX);
        Rz<<cos(phyZ),sin(phyZ),0,-sin(phyZ),cos(phyZ),0,0,0,1;
        position = Rz*Rx*xyzGK;
    } else {
        position = transfer*Vector2d(xk,yk);
    }


}

UbloxSolver::UbloxSolver(){
    rxyzt = Vector4d::Zero();
    for(int i = 0;i<32;i++){
        GPSSVs[i].svId = i+1;
        GPSSVs[i].type = SvInfo::GPS;
    }
    for(int i = 0;i<37;i++){
        BeiDouSVs[i].svId = i+1;
        BeiDouSVs[i].type = SvInfo::BeiDou;
    }
}

UbloxSolver::~UbloxSolver(){}


bool UbloxSolver::ParseRawData(char *message) {

    char* playload = message + 6;

    char* temp = playload;
    rcvtow = *(double*)temp;
    temp = playload + 11;
    numMeas = *(u_int8_t*)temp;
    if(0==numMeas)return false;

    for(u_int8_t n = 0;n<numMeas;n++){
        SvInfo* svTemp;
        int n32 = n*32;
        uint8_t gnssId = *(uint8_t *)(playload+36+n32);
        uint8_t svId = *(uint8_t *)(playload+37+n32);
        switch (gnssId){
            case 0:
                svTemp = &(GPSSVs[svId-1]);
                if(!useGPS)svTemp->open = false;
                break;
            case 3:
                svTemp = &(BeiDouSVs[svId-1]);
                if(!useBeiDou)svTemp->open = false;
                break;
        }
        visibleSvs.push_back(svTemp);
        svTemp->prMes = *(double*)(playload+16+n32);
        svTemp->cpMes = *(double*)(playload+24+n32);
        svTemp->doMes = *(float *)(playload+32+n32);

    }
    printf("start solving rawdata numMesa=%d\n",numMeas);
    //将星历拷贝一份防止计算时被改.
    vector<SvInfo>().swap(SvsForCalcu);
    for(int i=0;i<numMeas;i++)
    {
        SvInfo sv = *visibleSvs[i];
        printf("svsfor visiable:%d-%02d:%d,%d,%d\n",sv.type,sv.svId,sv.page1OK,sv.page2OK,sv.page3OK);
        if(sv.pageOK&&sv.open){
            SvsForCalcu.push_back(sv);
            printf("svsfor calcu\n");
        }
    }
    if(SvsForCalcu.size()<5){
        printf("calcu:Not enough Svs.\n");
        return false;
    }

//    thread PositionThread(LaunchPositionThread,this);
//    PositionThread.detach();
    solvePosition();
    return true;
}

bool UbloxSolver::ParseBstSubFrame(char *message) {
    char* playload = message + 6;

    uint8_t gnssId = *(uint8_t*)(playload);
    uint8_t svId = *(uint8_t*)(playload + 1);
    uint8_t numWords = *(uint8_t*)(playload+4);

    printf("Update subframe;;gnssid:%d,svid:%d\n",gnssId,svId);

    char* tmp = playload+8;

    if(10==numWords){
        uint32_t dwrds[10];
        for(int i=0;i<10;i++)   dwrds[i] = *(uint32_t*)(tmp+4*i);
        switch (gnssId){
            case 0:
                DecodeGpsBroadcast(dwrds,&GPSSVs[svId-1]);
                break;
            case 3:
                if(svId>5){
                    DecodeBeiDouBroadcastD1(dwrds,&BeiDouSVs[svId-1]);
                }
                else{
                    BeiDouSVs[svId-1].isBeiDouGEO = true;
                    DecodeBeiDouBroadcastD2(dwrds,&BeiDouSVs[svId-1]);
                }
                break;
        }
    }

}

bool UbloxSolver::solvePosition() {
    if(isCalculating)
        return false;
    //todo:
    //选星？
    int N = SvsForCalcu.size();
    double omegat,tempx,tempy,tempz;
    rxyzOld = rxyzt;
    Vector4d dtxyzt = Vector4d::Ones();
    MatrixXd pc(N,1);
    MatrixXd b(N,1);

    for(int i = 0;i< N;i++){
        SvInfo* sv= &SvsForCalcu[i];
        sv->CalcuECEF(rcvtow);
        //todo: calcu I,T,dtu
        double pci = sv->prMes + Light_speed * sv->tsDelta - sv->I - sv->T;
        pc(i) = pci;
        sv->PrintInfo(1);
        cout<<sv->position.norm()<<endl;

    }

    int numCalcu = 0;
    while (dtxyzt.norm()>0.1){
        numCalcu++;
        MatrixXd G(N,4);
        for(int i = 0;i< N;i++){
            SvInfo* sv= &SvsForCalcu[i];
            Vector3d rxyz = rxyzt.head(3);
            double r = (sv->position-rxyz).norm();
            //自转修正,这个我需要考虑一下是不是在这里处理
            omegat = r/Light_speed*Omega_e;
            MatrixXd earthRotate(3,3);
            earthRotate<<cos(omegat),sin(omegat),0,-sin(omegat),cos(omegat),0,0,0,1;
            Vector3d svPositionEarthRotate = earthRotate*sv->position;
            cout<<"positon = \n"<<sv->position<<endl;
            cout<<"earthRotate = \n"<<earthRotate<<endl;
            cout<<"positionRotate = \n"<<svPositionEarthRotate<<endl;
            r = (svPositionEarthRotate-rxyz).norm();
            cout<<"r = "<<r<<endl;
//            tempx = sv->position(0)*cos(omegat)+sv->position(1)*sin(omegat);
//            tempy = -sv->position(0)*sin(omegat)+sv->position(1)*cos(omegat);
//            tempz = sv->position(2);

            b(i) = pc(i)-r-rxyzt(3);  //这里的r需要考虑一下
            G(i,0) = (rxyzt(0)-svPositionEarthRotate(0))/r;  //x
            G(i,1) = (rxyzt(1)-svPositionEarthRotate(1))/r;  //y
            G(i,2) = (rxyzt(2)-svPositionEarthRotate(2))/r;  //z
            G(i,3) = 1;
        }
        MatrixXd GT = G.transpose();
//        dtxyzt = ((GT*G).inverse())*GT*b;
        dtxyzt = GT*b;
        dtxyzt = (GT*G).inverse()*dtxyzt;
        cout<<"G = "<<G<<endl;
        cout<<"GT = "<<GT<<endl;
        cout<<"b = "<<b<<endl;
        cout<<"dt = "<<dtxyzt<<endl;
        rxyzt += dtxyzt;
        cout<<"++++++++++++xyzt\n"<<rxyzt<<endl;
        cout<<"xyzNorm\n"<<rxyzt.head(3).norm()<<endl;
        XYZ2LLA(rxyzt.head(3),LLA);
        cout<<"++++++++++++LLA\n"<<LLA<<endl;
        if(numCalcu>30)return false;

    }

    XYZ2LLA(rxyzt.head(3),LLA);
    cout<<"++++++++++++rxyzt\n"<<rxyzt<<endl;
    cout<<"++++++++++++LLA\n"<<LLA<<endl;
    return true;

}

bool UbloxSolver::CalcuSvTime() {

}

bool UbloxSolver::DecodeBeiDouBroadcastD1(uint32_t *dwrds, SvInfo *sv) {
    uint32_t bit32,bit32x2[2];
    uint64_t bit64;
    //word0:
    if(1810!=Read1Word(dwrds[0],11))
        return false;
    int page = Read1Word(dwrds[0],3,17);
    sv->SOW = Read2Word(dwrds,8,20,12,2);

    printf(" Frame BeidouD1 svid:%d,page:,%d",sv->svId,page);
    switch (page){
        case 1:
            sv->page1OK = true;
            sv->WN = Read1Word(dwrds[2],13);
            sv->ino.a0 = ((int32_t) Read1Word(dwrds[4],8,8,true))*pow(2,-30);
            sv->ino.a1 = ((int32_t) Read1Word(dwrds[4],8,16,true))*pow(2,-27);
            sv->ino.a2 = ((int32_t) Read1Word(dwrds[5],8,2,true))*pow(2,-24);
            sv->ino.a3 = ((int32_t) Read1Word(dwrds[5],8,10,true))*pow(2,-24);
            sv->ino.b0 = ((int32_t) Read2Word(dwrds+5,6,18,2,2,true))*pow(2,11);
            sv->ino.b1 = ((int32_t) Read1Word(dwrds[6],8,4,true))*pow(2,14);
            sv->ino.b2 = ((int32_t) Read1Word(dwrds[6],8,12,true))*pow(2,16);
            sv->ino.b3 = ((int32_t) Read2Word(dwrds+6,4,20,4,2,true))*pow(2,16);

            sv->AODC = Read1Word(dwrds[1],5,15);

            sv->toc = Read2Word(dwrds+2,9,15,8)*8;
            sv->a0 = (int32_t)Read2Word(dwrds+7,7,17,17,2,true)*pow(2,-33);
            sv->a1 = (int32_t)Read2Word(dwrds+8,5,19,17,2,true)*pow(2,-50);
            sv->a2 = (int32_t)Read1Word(dwrds[7],11,6,true)*pow(2,-66);

            sv->TGD1 = (int32_t)Read1Word(dwrds[3],10,10,true)*0.1;
            sv->TGD2 = (int32_t)Read2Word(dwrds+3,4,20,6,2,true)*0.1;
            sv->AODE = Read1Word(dwrds[9],5,19);
            break;
        case 2:
            sv->page2OK = true;
            sv->orbit.toeF2 = Read1Word(dwrds[9],2,22)<<15<<3;
            sv->orbit.toe = sv->orbit.toeF2|sv->orbit.toeF3;
            sv->orbit.sq_a = Read2Word(dwrds+8,12,12,20)*pow(2,-19);
            sv->orbit.e = Read2Word(dwrds+4,10,14,22)*pow(2,-33);
            sv->orbit.dtn = (int32_t)Read2Word(dwrds+1,10,14,6,2,true)*pow(2,-43);
            sv->orbit.M0 = (int32_t)Read2Word(dwrds+3,20,4,12,2,true)*pow(2,-31);
            sv->orbit.Cuc = (int32_t)Read2Word(dwrds+2,16,8,2,2,true)*pow(2,-31);
            sv->orbit.Cus = (int32_t)Read1Word(dwrds[6],18,2,true)*pow(2,-31);
            sv->orbit.Crc = (int32_t)Read2Word(dwrds+6,4,20,14,2,true)*pow(2,-6);
            sv->orbit.Crs = (int32_t)Read2Word(dwrds+7,8,16,10,2,true)*pow(2,-6);
            break;
        case 3:
            sv->page3OK = true;
            sv->orbit.toeF3 = Read2Word(dwrds+1,10,14,5)<<3;
            sv->orbit.toe = sv->orbit.toeF2|sv->orbit.toeF3;
            sv->orbit.omega = (int32_t)Read2Word(dwrds+8,11,13,21,2,true)*pow(2,-31);
            sv->orbit.Omega0 = (int32_t)Read2Word(dwrds+7,21,3,11,2,true)*pow(2,-31);
            sv->orbit.OmegaDot = (int32_t)Read2Word(dwrds+4,11,13,13,2,true)*pow(2,-43);
            sv->orbit.i0 = (int32_t)Read2Word(dwrds+2,17,7,15,2,true)*pow(2,-31);
            sv->orbit.IDOT = (int32_t)Read2Word(dwrds+6,13,11,1,2,true)*pow(2,-43);
            sv->orbit.Cic = (int32_t)Read2Word(dwrds+3,7,17,11,2,true)*pow(2,-31);
            sv->orbit.Cis = (int32_t)Read2Word(dwrds+5,9,15,9,2,true)*pow(2,-31);
            break;
    }
    sv->pageOK = sv->page1OK && sv->page2OK && sv->page3OK;
}

bool UbloxSolver::DecodeBeiDouBroadcastD2(uint32_t *dwrds, SvInfo *sv) {
    printf(" Frame BeidouD2 page:");

}

bool UbloxSolver::DecodeGpsBroadcast(uint32_t *dwrds, SvInfo *sv) {
    //todo:
    printf(" Frame GPS  page:");

}

uint32_t UbloxSolver::Read1Word(uint32_t word, int length, int head, bool isInt) {
    if(isInt)
        return uint32_t (((int32_t)word)<<head>>(32-length));
    else
        return word<<head>>(32-length);
}

uint32_t UbloxSolver::Read2Word(uint32_t *word, int length0, int head0, int length1, int head1, bool isInt) {
    uint32_t high,low;
    if(isInt)
        high = uint32_t (((int32_t)word[0])<<head0>>(32-length0)<<length1);
    else
        high = word[0]<<head0>>(32-length0)<<length1;
    low = word[1]<<head1>>(32-length1);
    return  high|low;
}

//xyz坐标系转换成纬经高坐标系  input:定位后的XYZ   output：转换后的纬经高
bool UbloxSolver::XYZ2LLA(Vector3d XYZ,Vector3d &LLA) {
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




    return true;
}

