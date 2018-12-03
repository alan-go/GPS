#include "Svs.h"
#include "NtripRTK.h"
#include "EphemSp3.h"

SV::SV(int id):svId(id),SatH1(1),I(0),T(0),isBeiDouGEO(false),elevationAngle(0),tsDelta(0),ephemSp3(nullptr),trackCount(0){
    memset(bstEphemOK,0,10 * sizeof(int8_t));
    kal = Kalman(30,29,0);
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
        sv->ion = &ionKlob;
    }
    sysAll.push_back(sBds);

    for(int i=0;i<Ngps;i++){
        SV* sv = new GpsSV(i+1);
        sGps->table.push_back(sv);
        sv->ion = &ionKlob;
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
                sprintf(tip,"LackFrame");
                return false;
            }
            if(SatH1){
                sprintf(tip,"noHeal");

                return false;
            }
            break;
        case 1:
            //todo
            if(ephemSp3== nullptr){
                sprintf(tip,"noSp3Data");
                return false;
            }
            timeperiod = ephemSp3->dt * 5;
            if((time-ephemSp3->timeHead<timeperiod)||(ephemSp3->timeEnd-time<timeperiod)){
                sprintf(tip,"sp3TimeError");
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

bool SV::AddMmeasure(Measure *mesr) {
    if(!measureDat.empty()){
        Measure* latest = measureDat.front();
        mesr->cycle = latest->cycle;
        mesr->cycleP = latest->cycleP;
        if(mesr->time-latest->time>1.5)trackCount=0;
        else trackCount++;
    }
    measureDat.push_front(mesr);
};

void SV::PrintInfo(int printType) {
    //1:position and tsDelta
    //2:bst ephemeris
    switch (printType){
        case 0:
            printf("sv:%d,%02d  %s", type, svId,tip);
            break;
        case 1:
//            cout<<"+++++++SvPosition:"<<type<<","<<svId<<endl;
//            cout<<position<<endl;
            printf("+++SvPosition:%d,%d === %10f,%10f,%10f\n",type,svId,xyz(0),xyz(1),xyz(2));
            cout<<"norm="<<xyz.norm()<<endl;
            printf("+++svLLA === %lf, %lf, %lf\n",lla(0)*180/GPS_PI,lla(1)*180/GPS_PI,lla(2));
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

void SV::FPrintInfo(int printType){
    if(NULL==fpLog){
    char name[32];
    sprintf(name,"../log/SV/%d_%02d.txt",type,svId);
    fpLog = fopen(name,"w");
    }
    Measure *ms = measureDat.front();
    Measure *ms1 = measureDat[1];
    double dt = ms->time.tow-(ms1->time.tow);
    switch(printType){
        case 0:

            break;
        case 1:
            fprintf(fpLog,"%.4f,\t%.4f,%.4f,\t",ms->time.tow,ms->prMes,ms->stdevPr);//1,2
            fprintf(fpLog,"%.4f,%.4f,\t",ms->cpMes,ms->stdevCp);//3,4
            fprintf(fpLog,"%.4f,%.4f,\t",ms->doMes,ms->stdevDo);//5,6
            fprintf(fpLog,"%.4f,%d,\t",ms->lockTime/1000,ms->trkStat);//7,8
            fprintf(fpLog,"%.4f,%.4f,%.4f\n",ms->prMes-ms->cpMes,(ms->prMes-ms1->prMes)/dt,(ms->cpMes-ms1->cpMes)/dt);//9,10,11
            break;
        case 2:
            break;
    }
}
bool SV::CalcuECEF(double ts0) {
    double A = orbit.sqrtA*orbit.sqrtA;
    double sqt1_e2 = sqrt(1-orbit.e*orbit.e);
    double n0 = sqrt(M_miu/(A*A*A));
    double tk = ts0 - orbit.toe;
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
    double cosEk=cos(Ek),sinEk=sin(Ek);
    double rsinvk = sqt1_e2*sinEk;
    double rcosvk = cosEk-orbit.e;
    double rvk = 1-cosEk*orbit.e;
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
    double xk_ = rk * cos(uk);
    double yk_ = rk * sin(uk);
    double Omega_e_ = isBeiDouGEO?0:Omega_e;
    double Omegak = orbit.Omega0 + (orbit.OmegaDot - Omega_e_)*tk - Omega_e*orbit.toe;
    MatrixXd transfer(3,2);
    double cosik=cos(ik),sinik=sin(ik),cosOmk=cos(Omegak),sinOmk=sin(Omegak);
    transfer<<cosOmk,-cosik*sinOmk,sinOmk,cosik*cosOmk,0,sinik;
    Vector3d xyzGK = transfer*Vector2d(xk_,yk_);

    //calcu SV rate(according to book page64)
    //step12
    double EkDot = n/(1-orbit.e*cosEk);
    double phyKDot = sqt1_e2*EkDot/(1-orbit.e*cosEk);

    double dtukDot = 2*phyKDot*(orbit.Cus*cos2phy-orbit.Cuc*sin2phy);
    double dtrkDot = 2*phyKDot*(orbit.Crs*cos2phy-orbit.Crc*sin2phy);
    double dtikDot = 2*phyKDot*(orbit.Cis*cos2phy-orbit.Cic*sin2phy);

    double ukDot = phyKDot + dtukDot;
    double rkDot = A*orbit.e*EkDot*sinEk + dtrkDot;
    double ikDot = orbit.IDOT + dtikDot;
    double OmegakDot = orbit.OmegaDot - Omega_e_;

    double xkDot_ = rkDot*cos(uk)-rk*ukDot*sin(uk);
    double ykDot_ = rkDot*sin(uk)+rk*ukDot*cos(uk);

    double xDot = -xyzGK(1)*OmegakDot-(ykDot_*cosik-xyzGK(2)*ikDot)*sinOmk+xkDot_*cosOmk;
    double yDot = xyzGK(0)*OmegakDot+(ykDot_*cosik-xyzGK(2)*ikDot)*cosOmk+xkDot_*sinOmk;
    double zDot = ykDot_*sinik+yk_*ikDot*cosik;
    Vector3d vxyzGK(xDot,yDot,zDot);
    if(isBeiDouGEO){
        Matrix3d Rz,Rx,Rz_,S;
        double phyX, phyZ;
        phyX = -5.0/180.0*GPS_PI;
        phyZ = Omega_e * tk;
        double cosphyX = cos(phyX),sinphyX=sin(phyX),cosphyZ=cos(phyZ),sinphyZ=sin(phyZ);
        Rx<<1,0,0,0,cosphyX,sinphyX,0,-sinphyX,cosphyX;
        Rz<<cosphyZ,sinphyZ,0,-sinphyZ,cosphyZ,0,0,0,1;
        Rz_<<-Omega_e*sinphyZ,Omega_e*cosphyZ,0,-Omega_e*cosphyZ,-Omega_e*sinphyZ,0,0,0,0;
        xyz = Rz*Rx*xyzGK;
        vxyz = Rz_*Rx*xyzGK+Rz*Rx*vxyzGK;

    } else{
        xyz=xyzGK;
        vxyz=vxyzGK;
    }

//todo why do tgd and relativity here?
//TGD and relativity fix.
    double dtRelativity = 2.0*sqrt(M_miu)/(Light_speed*Light_speed)*Earth_ee*orbit.sqrtA*sin(Ek);
//    printf("dtRelativity=%.10f\n",dtRelativity);
    tsDelta -= dtRelativity;
    tsDelta -= TGD1;
    tsReal= ts0-tsDelta;
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
            ion->a0 = ((int32_t) Read1Word(dwrds[4],8,8,true))*pow(2,-30);
            ion->a1 = ((int32_t) Read1Word(dwrds[4],8,16,true))*pow(2,-27)/GPS_PI;
            ion->a2 = ((int32_t) Read1Word(dwrds[5],8,2,true))*pow(2,-24)/GPS_PI2;
            ion->a3 = ((int32_t) Read1Word(dwrds[5],8,10,true))*pow(2,-24)/GPS_PI3;
            ion->b0 = ((int32_t) Read2Word(dwrds[5],6,18,dwrds[6],2,2,true))*pow(2,11);
            ion->b1 = ((int32_t) Read1Word(dwrds[6],8,4,true))*pow(2,14)/GPS_PI;
            ion->b2 = ((int32_t) Read1Word(dwrds[6],8,12,true))*pow(2,16)/GPS_PI2;
            ion->b3 = ((int32_t) Read2Word(dwrds[6],4,20,dwrds[7],4,2,true))*pow(2,16)/GPS_PI3;

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
            ion->a0 = ((int32_t) Read2Word(dwrds[1],6,18,dwrds[2],2,2,true))*pow(2,-30);
            ion->a1 = ((int32_t) Read1Word(dwrds[2],8,4,true))*pow(2,-27)/GPS_PI;
            ion->a2 = ((int32_t) Read1Word(dwrds[2],8,12,true))*pow(2,-24)/GPS_PI2;
            ion->a3 = ((int32_t) Read2Word(dwrds[2],4,20,dwrds[3],4,2,true))*pow(2,-24)/GPS_PI3;
            ion->b0 = ((int32_t) Read1Word(dwrds[3],8,6,true))*pow(2,11);
            ion->b1 = ((int32_t) Read1Word(dwrds[3],8,14,true))*pow(2,14)/GPS_PI;
            ion->b2 = ((int32_t) Read2Word(dwrds[3],2,22,dwrds[4],6,2,true))*pow(2,16)/GPS_PI2;
            ion->b3 = ((int32_t) Read1Word(dwrds[4],8,8,true))*pow(2,16)/GPS_PI3;
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

int SV::CalcuelEvationAzimuth(Vector3d pos, Vector3d poslla) {
    Vector3d dtenu,dtxyz = xyzR - pos;
    XYZ2ENU(dtxyz,poslla,dtenu);
    elevationAngle = asin(dtenu(2)/dtenu.norm());
    azimuthAngle = atan2(dtenu(0),dtenu(1));
    if(azimuthAngle<0)azimuthAngle+=2*GPS_PI;
}

int SV::CalcuInoshphere(double elev, double azim,Vector3d LLA,double time) {
    double temp0 = Earth_a/(Earth_a+375000);
    double temp1 = temp0*cos(elev);
//    printf("ai,bi %e,%ef,%ef,%ef,\t %e,%e,%e,%e,\n",ion->a0,ino.a1,ino.a2,ino.a3,ino.b0,ino.b1,ino.b2,ino.b3 );
    double psi = GPS_PI/2 - elev - asin(temp1);
    double phyM = asin(sin(LLA(0))*cos(psi) + cos(LLA(0))*sin(psi)*cos(azim));
    double lambdaM = LLA(1)+asin(sin(psi)*sin(azim)/cos(phyM));
    double t = time + lambdaM*43200/GPS_PI;
    t = fmod(t,86400);
    if(t<0)t+=86400;

    double phyMpi = phyM/GPS_PI;
    double A2 = ion->a0 + phyMpi * (ion->a1 + phyMpi * (ion->a2 + phyMpi * ion->a3));
    if(A2<0)A2 = 0;
    double A4 = ion->b0 + phyMpi * (ion->b1 + phyMpi * (ion->b2 + phyMpi * ion->b3));
    if(A4<72000)A4 = 72000;
    if(A4>172800)A4 = 172800;

    double Iz_ = 5e-9;
    if(abs(t - 50400) < A4/4)
        Iz_ += A2*cos(2*GPS_PI*(t-50400)/A4);

    I = Iz_ / sqrt(1 - temp1*temp1);
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

int SV::CorrectIT(Vector3d xyz,Vector3d lla,double time) {
    CalcuelEvationAzimuth(xyz,lla);
//    printf("%d,%d..elevation = %lf, azim = %lf\n",type,svId,elevationAngle,azimuthAngle);
    CalcuTroposhphere(elevationAngle,azimuthAngle);
    CalcuInoshphere(elevationAngle,azimuthAngle,lla,time);
}



bool SV::ElevGood() {
    if(0==elevationAngle)
        return 1;
    elevGood = elevationAngle>0.17;
    if(!elevGood){
        sprintf(tip,"elevTooLow:%.1f",elevationAngle);
    }
    return elevGood;
}

bool SV::MeasureGood() {
    Measure *temp = measureDat.front();
    measureGood = temp->prMes<45e6&&temp->prMes>15e6;
    if(!measureGood){
        sprintf(tip,"prError:%.1f",temp->prMes);
    }
//    if (trackCount<5) {
//        sprintf(tip,"trackCount %d<5", trackCount);
//        return 0;
//    }

    return measureGood;
}

double SV::InterpRtkData(double time, int sigInd) {
    int InterpLength = 5;

    auto ite = rtkData.begin();
    double distanse = 100;

    for (; ite < rtkData.end() - 5; ite++) {
        distanse = time - (*ite)->rtktime;
        if(distanse>0)break;
    }
    if(distanse>2)
    {
        sprintf(tip," %drtk,cloest=%.1f",rtkData.size(),distanse);
        return 0;
    }
//    printf("time: %lf,%lf\n", time,(*ite)->rtktime);

    vector<double> times, prMes, cpMes;
    for(int i=0;i<InterpLength;i++){
        MSM4data* data = *(ite+i);
        MSM4Cell* cell = &data->sigData[sigInd];
        if(cell->cpMes!=0){
            times.push_back(data->rtktime);
            prMes.push_back(cell->prMes);
            cpMes.push_back(cell->cpMes);
//            printf("t=%f,pr=%f,cp = %f\n",data->rtktime,cell->prMes,cell->cpMes);
        }

    }
    if(times.size()<InterpLength-2)
    {
        printf("interp time line size = %d, not ok\n",times.size());
        return 0;
    }
    prInterp[sigInd] = InterpLine(time,times,prMes);
    cpInterp[sigInd] = InterpLine(time,times,cpMes);
//    printf("Interp:t=%f,pr=%f,cp = %f\n",time,prInterp[sigInd],cpInterp[sigInd]);

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

double SV::InterpMeasere(int len, int power, int begin) {
    if(measureDat.size()<(len+begin))
        return -1;
    MatrixXd G(len,power+1),W(len,len);
    VectorXd Ai(power+1),b(len);
    W.fill(0);
    double t0 = measureDat[begin]->time.tow;
    double dt1 = t0-measureDat[begin+1]->time.tow;
    if(dt1>3){
//        measureDat[0]->cycle=0;
        return -1;
    }
    for(int i=begin;i<begin+len;i++){
        Measure* ms = measureDat[i];
        double dt = t0-ms->time.tow;
        for(int n=0;n<=power;n++)G(i,n)=pow(dt,n);
        W(i,i)=1/pow(ms->stdevDo,2);
        b(i)=ms->doMes;
    }
    MatrixXd WG = W*G;
//    MatrixXd WGT = WG.transpose();
    Ai = WG.colPivHouseholderQr().solve(W*b);
    //积分
    double dDopler =0;
    for(int n=0;n<=power;n++){
        dDopler+=Ai(n)*pow(dt1,n+1)/(n+1);
    }


    /////////////debug
    double lambda = GetFreq(type,1,1);
    Measure* ms0 = measureDat[begin];
    Measure* ms1 = measureDat[begin+1];
    double dcp = ms0->cpMes-ms1->cpMes;
    double dpr = ms0->prMes-ms1->prMes;
//    dDopler = (ms0->doMes+ms1->doMes)*dt1/2;
//    if(abs(dcp-dDopler)>500)  {
//        return -1;
//    }
    double Cycle = ((dcp-dDopler)/lambda);
    double dCycle = fmod(Cycle+1e6,1);
    if(dCycle>0.5)dCycle-=1;
    double temp = (Cycle-dCycle);
    ms0->cycle+=temp;
    double pr_cpC0 = ms0->prMes-ms0->cpMes;
    double pr_cpC_ = ms0->prMes-ms0->cpMes-(ms0->cycle)*lambda;
    double pr_cpC1 = ms0->prMes-ms0->cpMes+(ms0->cycle)*lambda;
    fprintf(fpLog,"%.4f,\t",ms0->time.tow);//0,
    fprintf(fpLog,"%f,%f,%f,%f,%f,    ",dpr,dcp,dDopler,Cycle,dCycle);//12345
    fprintf(fpLog,"%f,%f,%f,%f,   ",(dpr)/dt1,(dDopler)/dt1,ms0->stdevPr,ms0->stdevDo);//6789
    fprintf(fpLog,"%f,%f,%f\n",pr_cpC0,pr_cpC_,pr_cpC1);//10 11
}

double SV::SmoothPr(int len, int begin) {
    double lambda = GetFreq(type,1,1);
    Measure* ms0 = measureDat[begin];
    Measure* ms1 = measureDat[begin+1];
    double t0 = ms0->time.tow,t1 = ms1->time.tow;
    double weight0 = 1/pow(ms0->stdevDo,2),weight1 = 1/pow(ms1->stdevDo,2);
    double dt01 = t0-t1;
    if(dt01<1){ms0->trackTime=ms1->trackTime+dt01;}
    else { return -1;}

    double dDoppler = (weight0*ms0->doMes+weight1*ms1->doMes)/(weight0+weight1)*dt01;
    double dcp = ms0->cpMes-ms1->cpMes;
    double dDoCp = dcp-dDoppler;
    ms0->cycleSlip = round(dDoCp/lambda);
    ms0->cycleRes = dDoCp/lambda-ms0->cycleSlip;

    int N = 0;
    double cycleTemp=0,diffSum=0,diffMean=0;
    for(int i=begin;;i++,N++){
        Measure* ms = measureDat[i];
        if(0==ms->trackTime||ms->stdevDo>0.5)break;
        if(t0-ms->time.tow>len)break;
        cycleTemp+=ms->cycleSlip;
        double diff = ms->cpMes-ms->prMes-cycleTemp*lambda;
        diffMean = (diff+diffMean*N)/(N+1);
    }

    double stdSum=0;
    for(int i=begin;i<N;i++){
        Measure* ms = measureDat[i];
        cycleTemp+=ms->cycleSlip;
        double diff = ms->cpMes-ms->prMes-cycleTemp*lambda;
        double std = pow(diff-diffMean,2);
        stdSum+=std;
    }

}


double SV::SmoothKalman(int len, int begin) {
    Measure* ms0 = measureDat[begin];
    Measure* ms1 = measureDat[begin+1];
    double t0 = ms0->time.tow,t1 = ms1->time.tow;
    double dt01 = t0-t1;

    if(dt01>1.5){
        kal.state=0;
        return -1;
    }
    double weight0 = 1/pow(ms0->stdevDo,2),weight1 = 1/pow(ms1->stdevDo,2);
    if(dt01<1){ms0->trackTime=ms1->trackTime+dt01;}
    else { return -1;}
    ms0->dDoppler = (weight0*ms0->doMes+weight1*ms1->doMes)/(weight0+weight1)*dt01;

    kal.state++;
    int N = kal.N,M=kal.M;
    if(kal.state<kal.N)return -2;
    if(kal.state==kal.N){
        kal.Init();
        for(int ni=0;ni<N;ni++)kal.x(ni)=measureDat[ni]->prMes;
    }

    for(int ni=0;ni<N;ni++){
        Measure* msi = measureDat[ni];
        kal.y(ni)=msi->prMes;
        kal.Rmm(ni,ni) = pow(msi->stdevPr,2);
        if(N-1==ni)continue;
        kal.y(ni+N)=msi->dDoppler;
        kal.Rmm(ni+N,ni+N)=pow(msi->stdevDo,2);
    }

    kal.x.tail(N-1)=kal.x.head(N-1);
    kal.x(0)+=ms0->dDoppler;
    kal.Hmn.fill(0);
    kal.Hmn.block(0,0,N,N)=MatrixXd::Identity(N,N);
    kal.Hmn.block(N,0,N-1,N-1)+=MatrixXd::Identity(N-1,N-1);
    kal.Hmn.block(N,1,N-1,N-1)-=MatrixXd::Identity(N-1,N-1);
//    cout<<"Rn:"<<endl<<kal.Rmm<<endl;
//    cout<<"P:"<<endl<<kal.Pnn<<endl;

    fprintf(fpLog,"%.4f,\t %.10f",ms0->time.tow,ms0->prMes);//0,
//    fprintf(fpLog,"%.4f,\t %.10f",ms0->time.tow,kal.x(0));//0,

    kal.Rectify();
    fprintf(fpLog,",%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",kal.x(0),(kal.x(0)-kal.x(1))/dt01,ms0->dDoppler/dt01);//0,
}

double SV::SmoothKalman0() {
    Measure* ms0 = measureDat[0];
    ms0->prCor = ms0->prMes;
    ms0->stdPrCor = ms0->stdevPr;
    if(measureDat.size()<10)
        return -1;
    Measure* ms1 = measureDat[1];
    double t0 = ms0->time.tow,t1 = ms1->time.tow;
    double dt01 = t0-t1;

    if(dt01>1.5){
        kal.state=0;
        return -1;
    }
    double weight0 = 1/pow(ms0->stdevDo,2),weight1 = 1/pow(ms1->stdevDo,2);
    if(dt01<1){ms0->trackTime=ms1->trackTime+dt01;}
    else { return -1;}
    ms0->dDoppler = (weight0*ms0->doMes+weight1*ms1->doMes)/(weight0+weight1)*dt01;
    double dcp = ms0->cpMes-ms1->cpMes;

    double lambda = GetFreq(type,1,1);
    double dDoCp = dcp-ms0->dDoppler;
    ms0->cycleSlip = round(dDoCp/lambda);
    ms0->cycleRes = dDoCp/lambda-ms0->cycleSlip;

    kal.state++;
    int N = kal.N,M=kal.M;
    if(kal.state<kal.N)return -2;

    kal.Pnn.fill(0);
    for(int ni=0;ni<N;ni++){
        Measure* msi = measureDat[ni];
        kal.x(ni)=msi->prMes;
        kal.Pnn(ni,ni)=pow(msi->stdevPr,2);

        if(N-1==ni)continue;
        kal.y(ni)=msi->dDoppler;
        kal.Rmm(ni,ni)=pow(msi->stdevDo,2);
    }


    kal.Hmn.fill(0);
    kal.Hmn.block(0,0,N-1,N-1)+=MatrixXd::Identity(N-1,N-1);
    kal.Hmn.block(0,1,N-1,N-1)-=MatrixXd::Identity(N-1,N-1);
//    cout<<"Rn:"<<endl<<kal.Rmm<<endl;
//    cout<<"P:"<<endl<<kal.Pnn<<endl;

    fprintf(fpLog,"%.4f,\t %.10f",ms0->time.tow,ms0->prMes);//0,

    kal.Rectify();
    ms0->prCor = kal.x(0);
    ms0->stdPrCor = sqrt(kal.Pnn(0,0));
//    cout<<"Pafter:"<<endl<<kal.Pnn<<endl;
    fprintf(fpLog,",%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",kal.x(0),ms0->stdPrCor,ms0->prMes-kal.x(0));//0,
}
double SV::DetectCycleSlip() {
    if(measureDat.size()<2)
        return -1;
    Measure* ms0 = measureDat[0];
    Measure* ms1 = measureDat[1];
    double t0 = ms0->time.tow,t1 = ms1->time.tow;
    double dt01 = t0-t1;

    if(dt01>1.2){
        ms0->trackTime=0;
        return -1;
    } else{
        ms0->trackTime=ms1->trackTime+dt01;
    }
    double weight0 = 1/pow(ms0->stdevDo,2),weight1 = 1/pow(ms1->stdevDo,2);
    ms0->dDoppler = (weight0*ms0->doMes+weight1*ms1->doMes)/(weight0+weight1)*dt01;
    double dcp = ms0->cpMes-ms1->cpMes;

    double lambda = GetFreq(type,1,1);
    double dDoCp = dcp-ms0->dDoppler;
    ms0->cycleSlip = round(dDoCp/lambda);
    ms0->cycleRes = dDoCp/lambda-ms0->cycleSlip;
    //调参？
//    ms0->cycleSlipQ = pow(10*dt01*dt01/lambda,2);
    ms0->cycleSlipQ = pow(50/lambda*ms0->stdevDo*dt01*dt01,2);
}
