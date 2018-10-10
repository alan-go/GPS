#include "Svs.h"
#include "NtripRTK.h"
#include "EphemSp3.h"

SV::SV(int id):svId(id),SatH1(1),I(0),T(0),isBeiDouGEO(false),elevationAngle(0),tsDelta(0),ephemSp3(nullptr),trackCount(0){
    memset(bstEphemOK,0,10 * sizeof(int8_t));
}
SV::SV(){}
SV::~SV(){}
BeiDouSV::BeiDouSV(int id):SV(id){type = SYS_BDS;}
BeiDouSV::~BeiDouSV(){}
GpsSV::GpsSV(int id):SV(id){type = SYS_GPS;}
GpsSV::~GpsSV(){}


void SvAll::SetOpen(bool bds, bool gps) {
    GetSys(SYS_BDS)->OpenClose(bds);
    GetSys(SYS_GPS)->OpenClose(gps);
}
SvAll::SvAll(){}

void SvAll::InitAlloc() {
    SvSys* sBds = new SvSys(SYS_BDS);
    SvSys* sGps = new SvSys(SYS_GPS);
    for (int i = 0; i < Nbds; ++i) {
        SV* sv = new BeiDouSV(i+1);
        sBds->table.push_back(sv);
        sv->isBeiDouGEO = i<5?true:false;
    }
    sysAll.push_back(sBds);

    for(int i=0;i<Ngps;i++){
        sGps->table.push_back(new GpsSV(i+1));
    }
    sysAll.push_back(sGps);
}

//todo get usedSys
int SvAll::AddToUsed(SV *sv) {
    svUsedAll.push_back(sv);
    GetUsedSys(sv->type)->table.push_back(sv);
}

SvAll::~SvAll(){}


bool SV::IsEphemOK(int ephemType, GnssTime time) {
    int bstEphem;
    double timeperiod;
    switch (ephemType){
        case 0:
            bstEphem = bstEphemOK[0]*bstEphemOK[1]*bstEphemOK[2];
            if(isBeiDouGEO){
                for(int i = 3;i<10;i++){
                    bstEphem*=bstEphemOK[i];
                }
            }
            if(!bstEphem){
                printf("___Frame not enough\n");
                return false;
            }
            if(SatH1){
                printf("___not Healthy\n");

                return false;
            }
            break;
        case 1:
            //todo
            if(ephemSp3== nullptr){
                printf("_____no sp3 record \n");
                return false;
            }
            timeperiod = ephemSp3->dt * 5;
            if((time-ephemSp3->timeHead<timeperiod)||(ephemSp3->timeEnd-time<timeperiod)){
                printf("sp3 time  period not ok\n");
                return false;
            }
            break;
    }

    return true;
}

bool SV::CalcuTime(double tow) {
    tsv = tow - measureDat.front()->prMes/Light_speed;
    tsReal = tsv;
    //此时的钟差是没有考虑相对论效应和 TGD的
    for(int i = 0;i<5;i++){
        tsDelta = a0+a1*(tsReal-toc)+a2*(tsReal-toc)*(tsReal-toc);
        tsReal =tsv-tsDelta;
    }
    return true;
}


void SV::PrintInfo(int printType) {
    //1:position and tsDelta
    //2:bst ephemeris
    switch (printType){
        case 0:
            printf("type %d,id %d\n", type, svId);
            break;
        case 1:
//            cout<<"+++++++SvPosition:"<<type<<","<<svId<<endl;
//            cout<<position<<endl;
            printf("+++SvPosition:%d,%d === %10f,%10f,%10f\n",type,svId,position(0),position(1),position(2));
            cout<<"norm="<<position.norm()<<endl;
            printf("+++svLLA === %lf, %lf, %lf\n",sLLA(0)*180/GPS_PI,sLLA(1)*180/GPS_PI,sLLA(2));
            printf("+++TGD === %.3f\n",TGD1*1e9);

            cout<<"tsDelta === "<<tsDelta*Light_speed<<" ,a0"<<a0<<endl;
            break;
        case 2:
            printf("\n%d,%d\n",type,svId);
            printf("%d,\t%.10e,\t%.10e,\t%.10e,\t\n",SOW,a0,a1,a2);
            printf("%d,\t%.10e,\t%.10e,\t%.10e,\t\n",orbit.IODE,orbit.Crs,orbit.dtn,orbit.M0);
            printf("%.10e,\t%.10e,\t%.10e,\t%.10e,\t\n",orbit.Cuc,orbit.e,orbit.Cus,orbit.sqrtA);
            printf("%d,\t%.10e,\t%.10e,\t%.10e,\t\n",orbit.toe,orbit.Cic,orbit.Omega0,orbit.Cis);
            printf("%.10e,\t%.10e,\t%.10e,\t%.10e,\t\n",orbit.i0,orbit.Crc,orbit.omega,orbit.OmegaDot);
            printf("%.10e,\t%d,\t%d,\t%d,\t\n",orbit.IDOT,0,WN,0);
            printf("%d,\t%d,\t%.10e,\t%d,\t\n",0,SatH1,TGD1,IODC);
            break;
        default:
            break;
    };
}

bool SV::CalcuECEF(double tow) {
    double A = orbit.sqrtA*orbit.sqrtA;
    double n0 = sqrt(M_miu/(A*A*A));
    double tk = tow - orbit.toe;
    while(tk > 302400)tk-=604800;
    while(tk < -302400)tk+=604800;
    double n = n0 + orbit.dtn;
    double Mk = orbit.M0 + n*tk;
    while (Mk<0.0)Mk+=2.0*GPS_PI;
    while (Mk>2.0*GPS_PI)Mk-=2.0*GPS_PI;
    double Ek = Mk,EkOld = Ek-1;
    while(abs(Ek-EkOld)>1e-8){
        EkOld = Ek;
        Ek = EkOld-(EkOld-orbit.e*sin(EkOld)-Mk)/(1.0-orbit.e*cos(EkOld));
    }
    //todo:
    //Make Ek within (0-2pi)?????
//    cout<<"EK:"<<Ek<<endl;
    double rsinvk = (sqrt(1-orbit.e*orbit.e)*sin(Ek));
    double rcosvk = cos(Ek)-orbit.e;
    double rvk = 1-cos(Ek)*orbit.e;
    double sinvk = rsinvk/rvk;
    double cosvk = rcosvk/rvk;
//    double tanvk = rsinvk/rcosvk;
    double vk = atan2(rsinvk,rcosvk);
//    if(cosvk<0&&sinvk<0)vk-=GPS_PI;
//    if(cosvk<0&&sinvk>0)vk+=GPS_PI;

    double phyk = vk + orbit.omega;
    double sin2phy = sin(2.0*phyk),cos2phy = cos(2.0*phyk);
    double dtuk = orbit.Cus*sin2phy + orbit.Cuc*cos2phy;
    double dtrk = orbit.Crs*sin2phy + orbit.Crc*cos2phy;
    double dtik = orbit.Cis*sin2phy + orbit.Cic*cos2phy;
    double uk = phyk + dtuk;
//    printf("phyk=%10f,dtuk=%.10f,dtrk=%.10f,dtik=%.10f\nuk=%.10f",phyk,dtuk,dtrk,dtik,uk);
    double rk = A*rvk + dtrk;
//    printf("rk=%.10f\n",rk);
    double ik = orbit.i0 + orbit.IDOT*tk + dtik;
    double xk = rk * cos(uk);
    double yk = rk * sin(uk);
    double Omegak = orbit.Omega0 + (orbit.OmegaDot - (isBeiDouGEO?0:Omega_e))*tk - Omega_e*orbit.toe;
    MatrixXd transfer(3,2);
    transfer<<cos(Omegak),-cos(ik)*sin(Omegak),sin(Omegak),cos(ik)*cos(Omegak),0,sin(ik);

    if(isBeiDouGEO){
        Vector3d xyzGK = transfer*Vector2d(xk,yk);
        double phyX = -5.0/180.0*GPS_PI;
        double phyZ = Omega_e * tk;
        Matrix3d Rz,Rx;
        Rx<<1,0,0,0,cos(phyX),sin(phyX),0,-sin(phyX),cos(phyX);
        Rz<<cos(phyZ),sin(phyZ),0,-sin(phyZ),cos(phyZ),0,0,0,1;
        position = Rz*Rx*xyzGK;
    } else {
        position = transfer*Vector2d(xk,yk);
    }

//TGD and relativity fix.
    double dtRelativity = 2.0*sqrt(M_miu)/(Light_speed*Light_speed)*Earth_ee*orbit.sqrtA*sin(Ek);
//    printf("dtRelativity=%.10f\n",dtRelativity);
    tsDelta -= dtRelativity;
    tsDelta -= TGD1;
    tsReal-=tsDelta;
}

void SvAll::UpdateEphemeris(char *subFrame) {
    char* playload = subFrame + 6;
    uint8_t gnssId = *(uint8_t*)(playload);
    uint8_t svId = *(uint8_t*)(playload + 1);
    uint8_t numWords = *(uint8_t*)(playload+4);
//    printf("Update subframe;;gnssid:%d,svid:%d\n",gnssId,svId);
    char* tmp = playload+8;
    if(10==numWords){
        uint32_t dwrds[10];
        for(int i=0;i<10;i++) {
            dwrds[i] = *(uint32_t*)(tmp+4*i);
        }
        GetSv(SysType(gnssId),svId)->DecodeSubFrame(dwrds);
    }
}

uint32_t SV::Read1Word(uint32_t word, int length, int head, bool isInt) {
    if(isInt)
        return uint32_t (((int32_t)word)<<head>>(32-length));
    else
        return word<<head>>(32-length);
}


uint32_t SV::Read2Word(uint32_t word0, int length0, int head0,
                                uint32_t word1,int length1, int head1, bool isInt) {
    uint32_t high,low;
    if(isInt)
        high = uint32_t (((int32_t)word0)<<head0>>(32-length0)<<length1);
    else
        high = word0<<head0>>(32-length0)<<length1;
    low = word1<<head1>>(32-length1);
    return  high|low;
}

uint32_t SV::Read3Word(uint32_t word0, int length0, int head0, uint32_t word1, int length1, int head1, uint32_t word2,
                       int length2, int head2, bool isInt) {
    uint32_t high,low;
//    high = Read1Word(word0,length0,head0,isInt)<<length1<<length2;
//    low = Read2Word(word1,length1,head1,word2,length2,head2, false);
    high = Read2Word(word0,length0,head0,word1,length1,head1,isInt)<<length2;
    low = word2<<head2>>(32-length2);
    return  high|low;
}

int GpsSV::DecodeSubFrame(uint32_t *dwrds) {
    int gpsFrameHead = Read1Word(dwrds[0],8,2);
    if(139!=gpsFrameHead){
        printf("SYS_GPS frame Head matching failed. head = %d\n",gpsFrameHead);
        return false;
    }
//    sv->SatH1 = Read1Word(dwrds[1],1,19);
    uint32_t AS = Read1Word(dwrds[1],1,20);
    if(1==AS){
        printf("This SYS_GPS Satellite is working on A-S mode.\n");
    }
    int frame = Read1Word(dwrds[1],3,21);
    bstEphemOK[frame-1] = 1;
    printf(" Frame SYS_GPS  frame:%d",frame);

    uint32_t L2,PCodeState;
    switch(frame){
        case 1:
            WN = Read1Word(dwrds[2],10,2);
            L2 = Read1Word(dwrds[2],2,12);
            URAI = Read1Word(dwrds[2],4,14);
            SatH1 = Read1Word(dwrds[2],6,18);
            IODC = Read2Word(dwrds[3],2,24,dwrds[7],8,2);
            PCodeState = Read1Word(dwrds[3],1,2);
            TGD1 = (int32_t)Read1Word(dwrds[6],8,18,true)*pow(2,-31);
            toc = Read1Word(dwrds[7],16,10)*pow(2,4);
            a2 = (int32_t)Read1Word(dwrds[8],8,2,true)*pow(2,-55);
            a1 = (int32_t)Read1Word(dwrds[8],16,10,true)*pow(2,-43);
            a0 = (int32_t)Read1Word(dwrds[9],22,2,true)*pow(2,-31);
            break;
        case 2:
            //todo:IODE judge?
            orbit.IODE = Read1Word(dwrds[2],8,2);
            orbit.Crs = (int32_t)Read1Word(dwrds[2],16,10,true)*pow(2,-5);
            orbit.dtn = (int32_t)Read1Word(dwrds[3],16,2,true)*pow(2,-43)*GPS_PI;
            orbit.M0 = (int32_t)Read2Word(dwrds[3],8,18,dwrds[4],24,2,true)*pow(2,-31)*GPS_PI;
            orbit.Cuc = (int32_t)Read1Word(dwrds[5],16,2,true)*pow(2,-29);
            orbit.e = Read2Word(dwrds[5],8,18,dwrds[6],24,2)*pow(2,-33);
            orbit.Cus = (int32_t)Read1Word(dwrds[7],16,2,true)*pow(2,-29);
            orbit.sqrtA = Read2Word(dwrds[7],8,18,dwrds[8],24,2)*pow(2,-19);
            orbit.toe = Read1Word(dwrds[9],16,2)*pow(2,4);
            break;
        case 3:
            //IODE
            orbit.IODE = Read1Word(dwrds[9],8,2);
            orbit.Cic = (int32_t)Read1Word(dwrds[2],16,2,true)*pow(2,-29);
            orbit.Omega0 = (int32_t)Read2Word(dwrds[2],8,18,dwrds[3],24,2,true)*pow(2,-31)*GPS_PI;
            orbit.Cis = (int32_t)Read1Word(dwrds[4],16,2,true)*pow(2,-29);
            orbit.i0 = (int32_t)Read2Word(dwrds[4],8,18,dwrds[5],24,2,true)*pow(2,-31)*GPS_PI;
            orbit.Crc = (int32_t)Read1Word(dwrds[6],16,2,true)*pow(2,-5);
            orbit.omega = (int32_t)Read2Word(dwrds[6],8,18,dwrds[7],24,2,true)*pow(2,-31)*GPS_PI;
            orbit.OmegaDot = (int32_t)Read1Word(dwrds[8],24,2,true)*pow(2,-43)*GPS_PI;
            orbit.IDOT = (int32_t)Read1Word(dwrds[9],14,10,true)*pow(2,-43)*GPS_PI;
            break;
        default:
            break;
    }
    return 1;
}

int BeiDouSV::DecodeD1(uint32_t *dwrds) {
    if(1810!=Read1Word(dwrds[0],11,2))
        return false;
    int frame = Read1Word(dwrds[0],3,17);
    bstEphemOK[frame-1] = 1;
    SOW = Read2Word(dwrds[0],8,20,dwrds[1],12,2);
//    printf(" Frame BeidouD1 svid:%d,frame:,%d\n",svId,frame);
    switch (frame){
        case 1:
            SatH1 = Read1Word(dwrds[1],1,14);
            URAI = Read1Word(dwrds[1],4,20);
            if(URAI)printf("\n\n\nUARI not ok = %d\n\n\n",URAI);
            WN = Read1Word(dwrds[2],13,2);
            ino.a0 = ((int32_t) Read1Word(dwrds[4],8,8,true))*pow(2,-30);
            ino.a1 = ((int32_t) Read1Word(dwrds[4],8,16,true))*pow(2,-27);
            ino.a2 = ((int32_t) Read1Word(dwrds[5],8,2,true))*pow(2,-24);
            ino.a3 = ((int32_t) Read1Word(dwrds[5],8,10,true))*pow(2,-24);
            ino.b0 = ((int32_t) Read2Word(dwrds[5],6,18,dwrds[6],2,2,true))*pow(2,11);
            ino.b1 = ((int32_t) Read1Word(dwrds[6],8,4,true))*pow(2,14);
            ino.b2 = ((int32_t) Read1Word(dwrds[6],8,12,true))*pow(2,16);
            ino.b3 = ((int32_t) Read2Word(dwrds[6],4,20,dwrds[7],4,2,true))*pow(2,16);

            AODC = Read1Word(dwrds[1],5,15);

            toc = Read2Word(dwrds[2],9,15,dwrds[3],8,2)*8;
            a0 = (int32_t)Read2Word(dwrds[7],7,17,dwrds[8],17,2,true)*pow(2,-33);
            a1 = (int32_t)Read2Word(dwrds[8],5,19,dwrds[9],17,2,true)*pow(2,-50);
            a2 = (int32_t)Read1Word(dwrds[7],11,6,true)*pow(2,-66);

            TGD1 = (int32_t)Read1Word(dwrds[3],10,10,true)*1e-10;
            TGD2 = (int32_t)Read2Word(dwrds[3],4,20,dwrds[4],6,2,true)*1e-10;
            orbit.AODE = Read1Word(dwrds[9],5,19);
            break;
        case 2:
            orbit.toeHigh = Read1Word(dwrds[9],2,22)<<15;
            orbit.toe = (orbit.toeHigh|orbit.toeLow)*8;
            orbit.sqrtA = Read2Word(dwrds[8],12,12,dwrds[9],20,2)*pow(2,-19);
            orbit.e = Read2Word(dwrds[4],10,14,dwrds[5],22,2)*pow(2,-33);
            orbit.dtn = (int32_t)Read2Word(dwrds[1],10,14,dwrds[2],6,2,true)*pow(2,-43)*GPS_PI;
            orbit.M0 = (int32_t)Read2Word(dwrds[3],20,4,dwrds[4],12,2,true)*pow(2,-31)*GPS_PI;
            orbit.Cuc = (int32_t)Read2Word(dwrds[2],16,8,dwrds[3],2,2,true)*pow(2,-31);
            orbit.Cus = (int32_t)Read1Word(dwrds[6],18,2,true)*pow(2,-31);
            orbit.Crc = (int32_t)Read2Word(dwrds[6],4,20,dwrds[7],14,2,true)*pow(2,-6);
            orbit.Crs = (int32_t)Read2Word(dwrds[7],8,16,dwrds[8],10,2,true)*pow(2,-6);
            break;
        case 3:
            orbit.toeLow = Read2Word(dwrds[1],10,14,dwrds[2],5,2);
            orbit.toe = (orbit.toeHigh|orbit.toeLow)*8;
            orbit.omega = (int32_t)Read2Word(dwrds[8],11,13,dwrds[9],21,2,true)*pow(2,-31)*GPS_PI;
            orbit.Omega0 = (int32_t)Read2Word(dwrds[7],21,3,dwrds[8],11,2,true)*pow(2,-31)*GPS_PI;
            orbit.OmegaDot = (int32_t)Read2Word(dwrds[4],11,13,dwrds[5],13,2,true)*pow(2,-43)*GPS_PI;
            orbit.i0 = (int32_t)Read2Word(dwrds[2],17,7,dwrds[3],15,2,true)*pow(2,-31)*GPS_PI;
            orbit.IDOT = (int32_t)Read2Word(dwrds[6],13,11,dwrds[7],1,2,true)*pow(2,-43)*GPS_PI;
            orbit.Cic = (int32_t)Read2Word(dwrds[3],7,17,dwrds[4],11,2,true)*pow(2,-31);
            orbit.Cis = (int32_t)Read2Word(dwrds[5],9,15,dwrds[6],9,2,true)*pow(2,-31);
            break;
        default:
            break;
    }
}

int BeiDouSV::DecodeD2Frame1(uint32_t *dwrds) {
    //todo
//    return 0;
    if(1810!=Read1Word(dwrds[0],11,2))return -1;
    int frame = Read1Word(dwrds[0],3,17);
    if(1!=frame)return -1;
    SOW = Read2Word(dwrds[0],8,20,dwrds[1],12,2);
    int Pnum1 = Read1Word(dwrds[1],4,14);
    bstEphemOK[Pnum1-1] = 1;

//    printf(" Frame BeidouD2 svid:%d,frame1,page:,%d",svId,Pnum1);
    switch(Pnum1){
        case 1:
            SatH1 = Read1Word(dwrds[1],1,18);
            AODC = Read1Word(dwrds[1],5,19);
            URAI = Read1Word(dwrds[2],4,2);
            if(URAI)printf("\n\n\nUARI not ok = %d\n\n\n",URAI);
            WN = Read1Word(dwrds[2],13,6);
            toc = Read2Word(dwrds[2],5,19,dwrds[3],12,2)*8;
            TGD1 = (int32_t)Read1Word(dwrds[3],10,14,true)*1e-10;
            TGD2 = (int32_t)Read1Word(dwrds[4],10,2,true)*1e-10;
            break;
        case 2:
            ino.a0 = ((int32_t) Read2Word(dwrds[1],6,18,dwrds[2],2,2,true))*pow(2,-30);
            ino.a1 = ((int32_t) Read1Word(dwrds[2],8,4,true))*pow(2,-27);
            ino.a2 = ((int32_t) Read1Word(dwrds[2],8,12,true))*pow(2,-24);
            ino.a3 = ((int32_t) Read2Word(dwrds[2],4,20,dwrds[3],4,2,true))*pow(2,-24);
            ino.b0 = ((int32_t) Read1Word(dwrds[3],8,6,true))*pow(2,11);
            ino.b1 = ((int32_t) Read1Word(dwrds[3],8,14,true))*pow(2,14);
            ino.b2 = ((int32_t) Read2Word(dwrds[3],2,22,dwrds[4],6,2,true))*pow(2,16);
            ino.b3 = ((int32_t) Read1Word(dwrds[4],8,8,true))*pow(2,16);
            break;
        case 3:
            a0 = (int32_t)Read2Word(dwrds[3],12,12,dwrds[4],12,2,true)*pow(2,-33);
            a1High = Read1Word(dwrds[4],4,14,true)<<18;
            a1 = (int32_t)(a1High|a1Low)*pow(2,-50);
            break;
        case 4:
            a1Low = Read2Word(dwrds[1],6,18,dwrds[2],12,2);
            a1 = (int32_t)(a1High|a1Low)*pow(2,-50);
            a2 = (int32_t)Read2Word(dwrds[2],10,14,dwrds[3],1,2,true)*pow(2,-66);
            orbit.AODE = Read1Word(dwrds[3],5,3);
            orbit.dtn = (int32_t)Read1Word(dwrds[3],16,8,true)*pow(2,-43)*GPS_PI;
            orbit.CucHigh = Read1Word(dwrds[4],14,2,true)<<4;
            orbit.Cuc = (int32_t)(orbit.CucHigh|orbit.CucLow)*pow(2,-31);
            break;
        case 5:
            orbit.CucLow = Read1Word(dwrds[1],4,18);
            orbit.Cuc = (int32_t)(orbit.CucHigh|orbit.CucLow)*pow(2,-31);
            orbit.M0 = (int32_t)(Read3Word(dwrds[1],2,22,dwrds[2],22,2,dwrds[3],8,2,true))*pow(2,-31)*GPS_PI;
            orbit.Cus = (int32_t)Read2Word(dwrds[3],14,10,dwrds[4],4,2,true)*pow(2,-31);
            orbit.eHigh = Read1Word(dwrds[4],10,6)<<22;
            orbit.e = (orbit.eHigh|orbit.eLow)*pow(2,-33);
            break;
        case 6:
            orbit.eLow = Read2Word(dwrds[1],6,18,dwrds[2],16,2);
            orbit.e = (orbit.eHigh|orbit.eLow)*pow(2,-33);
            orbit.sqrtA = Read3Word(dwrds[2],6,18,dwrds[3],22,2,dwrds[4],4,2)*pow(2,-19);
            orbit.CicHigh = Read1Word(dwrds[4],10,6,true)<<8;
            orbit.Cic = (int32_t)(orbit.CicHigh|orbit.CicLow)*pow(2,-31);
            break;
        case 7:
            orbit.CicLow = Read2Word(dwrds[1],6,18,dwrds[2],2,2);
            orbit.Cic = (int32_t)(orbit.CicHigh|orbit.CicLow)*pow(2,-31);
            orbit.Cis = (int32_t)Read1Word(dwrds[2],18,4,true)*pow(2,-31);
            orbit.toe = Read2Word(dwrds[2],2,22,dwrds[3],15,2)*8;
            orbit.i0High = Read2Word(dwrds[3],7,17,dwrds[4],14,2,true)<<11;
            orbit.i0 = (int32_t)(orbit.i0High|orbit.i0Low)*pow(2,-31)*GPS_PI;
            break;
        case 8:
            orbit.i0Low = Read2Word(dwrds[1],6,18,dwrds[2],5,2);
            orbit.i0 = (int32_t)(orbit.i0High|orbit.i0Low)*pow(2,-31)*GPS_PI;
            orbit.Crc = (int32_t)Read2Word(dwrds[2],17,7,dwrds[3],1,2,true)*pow(2,-6);
            orbit.Crs = (int32_t)Read1Word(dwrds[3],18,3,true)*pow(2,-6);
            orbit.OmegaDotHigh = Read2Word(dwrds[3],3,21,dwrds[4],16,2,true)<<5;
            orbit.OmegaDot = (int32_t)(orbit.OmegaDotHigh|orbit.OmegaDotLow)*pow(2,-43)*GPS_PI;
            break;
        case 9:
            orbit.OmegaDotLow = Read1Word(dwrds[1],5,18);
            orbit.OmegaDot = (int32_t)(orbit.OmegaDotHigh|orbit.OmegaDotLow)*pow(2,-43)*GPS_PI;
            orbit.Omega0 = (int32_t)Read3Word(dwrds[1],1,23,dwrds[2],22,2,dwrds[3],9,2,true)*pow(2,-31)*GPS_PI;
            orbit.omegaHigh = Read2Word(dwrds[3],13,11,dwrds[4],14,2,true)<<5;
            orbit.omega = (int32_t)(orbit.omegaHigh|orbit.omegaLow)*pow(2,-31)*GPS_PI;
            break;
        case 10:
            orbit.omegaLow = Read1Word(dwrds[1],5,18);
            orbit.omega = (int32_t)(orbit.omegaHigh|orbit.omegaLow)*pow(2,-31)*GPS_PI;
            orbit.IDOT = (int32_t)Read2Word(dwrds[1],1,23,dwrds[2],13,2,true)*pow(2,-43)*GPS_PI;
            break;
        default:
            break;
    }
    return 1;
}

int SV::CalcuelEvationAzimuth(Vector3d receiverPosition, Vector3d LLA) {
    double x = receiverPosition(0);
    double y = receiverPosition(1);
    double z = receiverPosition(2);
    double r = sqrt(x*x+y*y);
    double sinB = y/r;
    double cosB = x/r;
    double sinA = sin(LLA(1));
    double cosA = cos(LLA(1));
    Matrix<double ,3,3>SS;
    SS<<-sinB,cosB,0,-sinA*cosB,-sinA*sinB,cosA,cosA*cosB,cosA*sinB,sinA;

    Vector3d dtxyz = position - receiverPosition;
    Vector3d dtenu = SS * dtxyz;
    elevationAngle = asin(dtenu(2)/dtenu.norm());
    azimuthAngle = atan(dtenu(0)/dtenu(1));
}

int SV::CalcuInoshphere(double elev, double azim,Vector3d LLA,double time) {
    double temp = Earth_a/(Earth_a+375000);
    double phy = GPS_PI/2 - elev - asin(temp*cos(elev));
    double phyM = asin(sin(LLA(1))*cos(phy) + cos(LLA(0))*sin(phy)*cos(azim));
    double lambdaM = LLA(0)+asin(sin(phy)*sin(azim)/cos(phyM));
    double t = time + lambdaM*43200/GPS_PI;
    t = fmod(t,86400);
    if(t<0)t+=86400;

    double phyMpi = phyM/GPS_PI;
    double A2 = ino.a0 + phyMpi * (ino.a1 + phyMpi * (ino.a2 + phyMpi * ino.a3));
    if(A2<0)A2 = 0;
    double A4 = ino.b0 + phyMpi * (ino.b1 + phyMpi * (ino.b2 + phyMpi * ino.b3));
    if(A4<72000)A4 = 72000;
    if(A4>172800)A4 = 172800;

    double Iz_ = 5e-9;
    if(abs(t - 50400) < A4/4)
        Iz_ += A2*cos(2*GPS_PI*(t-50400)/A4);

    I = Iz_ / sqrt(1 - temp*temp);
    I *= Light_speed;

    return 0;
}

int SV::CalcuTroposhphere(double elev, double azim) {
    if(elev<0.1116)
        T=20;
    else
        T=2.47/(sin(elev)+0.0121);
    return 0;
}

int SV::CorrectIT(Vector3d pos,double time) {
    Vector3d LLA;
    XYZ2LLA(pos,LLA);
    CalcuelEvationAzimuth(pos,LLA);
    printf("elevation = %lf, azim = %lf\n",elevationAngle,azimuthAngle);
    CalcuTroposhphere(elevationAngle,azimuthAngle);
    CalcuInoshphere(elevationAngle,azimuthAngle,LLA,time);
}



bool SV::ElevGood() {
    if(0==elevationAngle)
        return 1;
    elevGood = elevationAngle>0.17;
    if(!elevGood){
        printf("\nelev angle bad sv:%d,%02d,angle:%lf\n",type,svId,elevationAngle);
    }
    return elevGood;
}

bool SV::MeasureGood() {
    Measure *temp = measureDat.front();
    measureGood = temp->prMes<45e6&&temp->prMes>15e6;
    if(!measureGood){
        printf("\npr measure bad sv:%d,%02d,angle:%lf\n",type,svId,temp->prMes);
    }
    if (trackCount<5) {
        printf("trackingTime %d < 5\n", trackCount);
        return 0;
    }

    return measureGood;
}

double SV::InterpRtkData(double time, int sigInd) {
    int InterpLength = 5;

    auto ite = rtkData.begin();
    double distanse = 100;
//    for (; ite<rtkData.end()-5 ; ite++) {
//        int temp = abs((*ite)->rtktime-time);
//        if(distanse>temp){
//            distanse = temp;
//        }
//        if(distanse<0.5)break;
//    }
//    if(distanse>=0.5)
//    {
//        printf("rtk data num = %d, not enough,distanse=%lf\n",rtkData.size(),distanse);
//        return 0;
//    }

    for (; ite < rtkData.end() - 5; ite++) {
        distanse = time - (*ite)->rtktime;
        if(distanse>0)break;
    }
    if(distanse>2)
    {
        printf("rtk data num size = %d, not enough,closeest distanse=%lf\n",rtkData.size(),distanse);
        return 0;
    }
    printf("time: %lf,%lf\n", time,(*ite)->rtktime);

    vector<double> times, prMes, cpMes;
    for(int i=0;i<InterpLength;i++){
        MSM4data* data = *(ite+i);
        MSM4Cell* cell = &data->sigData[sigInd];
        if(cell->cpMes!=0){
            times.push_back(data->rtktime);
            prMes.push_back(cell->prMes);
            cpMes.push_back(cell->cpMes);
            printf("t=%f,pr=%f,cp = %f\n",data->rtktime,cell->prMes,cell->cpMes);
        }

    }
    if(times.size()<InterpLength-2)
    {
        printf("interp time line size = %d, not ok\n",times.size());
        return 0;
    }
    prInterp[sigInd] = InterpLine(time,times,prMes);
    cpInterp[sigInd] = InterpLine(time,times,cpMes);
    printf("Interp:t=%f,pr=%f,cp = %f\n",time,prInterp[sigInd],cpInterp[sigInd]);

    return prInterp[sigInd];
}

double SV::InterpLine(double xi, vector<double> &x, vector<double> &y) {
    if(x.size()!=y.size()) return 0;
    int N = x.size();
    MatrixXd x2(N,2);
    MatrixXd y1(N,1);
    for(int i=0;i<N;i++){
        x2(i,0) = 1;
        x2(i,1) = x[i]-x[0];
        y1(i) = y[i]-y[0];
//        printf("xi,yi = %f,%f\n",x2(i,1),y1(i));
    }
    MatrixXd x2T = x2.transpose();
    Vector2d ab = (x2T*x2).inverse()*(x2T*y1);
//    printf("ab = %f,%f\n",ab(0),ab(1));
    return ab(0) + ab(1)*(xi-x[0]) + y[0];
}

bool SV::IsMaskOn() {
    switch (type){
        case SYS_BDS:
            break;
    }
}
