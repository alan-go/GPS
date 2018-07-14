#include "SVs.h"

SV::SV(){
    page1OK = page2OK = page3OK = false;
    I = T = 0;
}
SV::~SV(){}

SVs::SVs(GNSS* gnss):gnss(gnss){
    for(int i = 0; i < 5; i++){
        svBeiDous[i].isBeiDouGEO = true;
    }
}

bool SV::JudgeUsable(bool useBeiDou, bool useGps) {
    switch (type){
        case GPS:
            if(!useGps)return false;
            break;
        case BeiDou:
            if(!useBeiDou)return false;
            break;
    }
    if(!(pageOK&&SatH1))return false;
    return true;
}

bool SV::CalcuTime(double rcvtow) {
    tsv = rcvtow - prMes/Light_speed;
    tsReal = tsv;
    //此时的钟差是没有考虑相对论效应和 TGD的
    for(int i = 0;i<5;i++){
        tsDelta = a0+a1*(tsReal-toc)+a2*(tsReal-toc)*(tsReal-toc);
        tsReal =tsv-tsDelta;
    }
    return true;
}


void SV::PrintInfo(int printType) {
    if(1==printType){
        cout<<"+++++++SvPosition:"<<type<<","<<svId<<endl;
        cout<<position<<endl;
        cout<<"norm="<<position.norm()<<endl;
        cout<<"tsDelta"<<tsDelta<<"a0"<<a0<<endl;
    }
}

bool SV::CalcuECEF(double rcvtow) {
    CalcuTime(rcvtow);
    double A = orbit.sqrtA*orbit.sqrtA;
    double n0 = sqrt(M_miu/(A*A*A));
    double tk = tsReal - orbit.toe;
    if(tk > 302400)tk-=604800;
    if(tk < -302400)tk+=604800;
    double n = n0 + orbit.dtn;
    double Mk = orbit.M0 + n*tk;
    while (Mk<0)Mk+=2*GPS_PI;
    while (Mk>2*GPS_PI)Mk-=2*GPS_PI;
    double Ek = Mk,EkOld = Ek-1;
    while(abs(Ek-EkOld)>1e-8){
        EkOld = Ek;
        Ek = EkOld-(EkOld-orbit.e*sin(EkOld)-Mk)/(1-orbit.e*cos(EkOld));
    }
    //todo:
    //Make Ek within (0-2pi)?????
    cout<<"EK:"<<Ek<<endl;
    double rsinvk = (sqrt(1-orbit.e*orbit.e)*sin(Ek));
    double rcosvk = cos(Ek)-orbit.e;
    double rvk = 1-cos(Ek)*orbit.e;
    double sinvk = rsinvk/rvk;
    double cosvk = rcosvk/rvk;
    double tanvk = rsinvk/rcosvk;
    double vk = atan(tanvk);
    if(cosvk<0&&sinvk<0)vk-=GPS_PI;
    if(cosvk<0&&sinvk>0)vk+=GPS_PI;
    cout<<"vk="<<vk<<endl;
    printf("vk = %10f\n",vk);
    double phyk = vk + orbit.omega;
    double sin2phy = sin(2*phyk),cos2phy = cos(2*phyk);
    double dtuk = orbit.Cus*sin2phy + orbit.Cuc*cos2phy;
    double dtrk = orbit.Crs*sin2phy + orbit.Crc*cos2phy;
    double dtik = orbit.Cis*sin2phy + orbit.Cic*cos2phy;
    double uk = phyk + dtuk;
    printf("phyk=%10f,dtuk=%.10f,dtrk=%.10f,dtik=%.10f\nuk=%.10f",phyk,dtuk,dtrk,dtik,uk);
    double rk = A*(1-orbit.e*cos(Ek)) + dtrk;
    printf("rk=%.10f\n",rk);
    double ik = orbit.i0 + orbit.IDOT*tk + dtik;
    double xk = rk * cos(uk);
    double yk = rk * sin(uk);
    double Omegak = orbit.Omega0 + (orbit.OmegaDot - (isBeiDouGEO?0:Omega_e)) * tk - Omega_e * orbit.toe;
    cout<<"OmegaK="<<Omegak<<endl;
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
//NO TGD here.
    double dtRelativity = -2*sqrt(M_miu)/(Light_speed*Light_speed)*orbit.sqrtA*sin(Ek);
    printf("dtRelativity=%.10f\n",dtRelativity);
    tsDelta += dtRelativity;
}

void SVs::UpdateEphemeris(char *subFrame) {
    char* playload = subFrame + 6;

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
                svGpss[svId-1].DecodeSubFrame(dwrds);
                break;
            case 3:
                svBeiDous[svId-1].DecodeSubFrame(dwrds);
                break;
        }
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

int GpsSV::DecodeSubFrame(uint32_t *dwrds) {
    printf(" Frame GPS  page:");
    int gpsFrameHead = Read1Word(dwrds[0],8,2);
    if(139!=gpsFrameHead){
        printf("GPS frame Head matching failed. head = %d\n",gpsFrameHead);
        return false;
    }
//    sv->SatH1 = Read1Word(dwrds[1],1,19);
    uint32_t SA = Read1Word(dwrds[1],1,20);
    if(1==SA){
        printf("This GPS Satellite is working on SA mode.\n");
    }
    int page = Read1Word(dwrds[1],3,21);
    uint32_t L2,PCodeState;
    switch(page){
        case 1:
            page1OK = true;
            WN = Read1Word(dwrds[2],10,2);
            L2 = Read1Word(dwrds[2],2,12);
            URAI = Read1Word(dwrds[2],4,14);
            SatH1 = Read1Word(dwrds[2],6,18);
            IODC = Read2Word(dwrds[3],2,24,dwrds[7],8,2);
            PCodeState = Read1Word(dwrds[3],1,2);
            TGD = (int32_t)Read1Word(dwrds[6],8,18,true)*pow(2,-31);
            toc = Read1Word(dwrds[7],16,10)*pow(2,4);
            a2 = (int32_t)Read1Word(dwrds[8],8,2,true)*pow(2,-55);
            a1 = (int32_t)Read1Word(dwrds[8],16,10,true)*pow(2,-43);
            a0 = (int32_t)Read1Word(dwrds[9],22,2,true)*pow(2,-31);
            break;
        case 2:
            page2OK = true;
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
            page3OK = true;
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
    pageOK = page1OK && page2OK && page3OK;
    return 1;
}

int BeiDouSV::DecodeD1(uint32_t *dwrds) {
    if(1810!=Read1Word(dwrds[0],11,2))
        return false;
    int page = Read1Word(dwrds[0],3,17);
    SOW = Read2Word(dwrds[0],8,20,dwrds[1],12,2);

    printf(" Frame BeidouD1 svid:%d,page:,%d",svId,page);
    switch (page){
        case 1:
            page1OK = true;
            SatH1 = Read1Word(dwrds[1],1,14);
            URAI = Read1Word(dwrds[1],4,20);
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

            TGD1 = (int32_t)Read1Word(dwrds[3],10,10,true)*0.1;
            TGD2 = (int32_t)Read2Word(dwrds[3],4,20,dwrds[4],6,2,true)*0.1;
            orbit.AODE = Read1Word(dwrds[9],5,19);
            break;
        case 2:
            ///
            page2OK = true;
            orbit.toeF2 = Read1Word(dwrds[9],2,22)<<15<<3;
            orbit.toe = orbit.toeF2|orbit.toeF3;
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
            page3OK = true;
            orbit.toeF3 = Read2Word(dwrds[1],10,14,dwrds[2],5,2)<<3;
            orbit.toe = orbit.toeF2|orbit.toeF3;
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
    pageOK = page1OK && page2OK && page3OK;
}

int BeiDouSV::DecodeD2(uint32_t *dwrds) {
    //todo
}