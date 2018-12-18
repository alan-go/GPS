//
// Created by alan on 18-12-12.
//

#include "EphemBst.h"

EphemBst::EphemBst() {}

EphemBst::EphemBst(SV *sv_):sv(sv_) {
//    memset(framOk,0,10 * sizeof(int8_t));
}

int EphemBst::CalcuTs(GnssTime ts0) {
    sv->ts0=sv->ts= ts0;
    //此时的钟差是没有考虑相对论效应和 TGD的
    for(int i = 0;i<5;i++){
        double dt = sv->ts.tow-clk.toc;
        sv->tsdt= clk.a0+clk.a1*dt+clk.a2*dt*dt;
        sv->ts =ts0-sv->tsdt;
    }
    return 0;
}

int EphemBst::CalcuECEF(GnssTime time) {
    double A = orbit.sqrtA*orbit.sqrtA;
    double sqt1_e2 = sqrt(1-orbit.e*orbit.e);
    double n0 = sqrt(M_miu/(A*A*A));
    double tk = time.tow - orbit.toe;
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
    double Omega_e_ = sv->isBeiDouGEO?0:Omega_e;
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
    if(sv->isBeiDouGEO){
        Matrix3d Rz,Rx,Rz_,S;
        double phyX, phyZ;
        phyX = -5.0/180.0*GPS_PI;
        phyZ = Omega_e * tk;
        double cosphyX = cos(phyX),sinphyX=sin(phyX),cosphyZ=cos(phyZ),sinphyZ=sin(phyZ);
        Rx<<1,0,0,0,cosphyX,sinphyX,0,-sinphyX,cosphyX;
        Rz<<cosphyZ,sinphyZ,0,-sinphyZ,cosphyZ,0,0,0,1;
        Rz_<<-Omega_e*sinphyZ,Omega_e*cosphyZ,0,-Omega_e*cosphyZ,-Omega_e*sinphyZ,0,0,0,0;
        sv->xyz = Rz*Rx*xyzGK;
        sv->vxyz = Rz_*Rx*xyzGK+Rz*Rx*vxyzGK;

    } else{
        sv->xyz=xyzGK;
        sv->vxyz=vxyzGK;
    }

//todo why do tgd and relativity here?
//TGD and relativity fix.
    double dtRelativity = 2.0*sqrt(M_miu)/(Light_speed*Light_speed)*Earth_ee*orbit.sqrtA*sin(Ek);
//    printf("dtRelativity=%.10f\n",dtRelativity);
    sv->tsdt-= dtRelativity;
    sv->tsdt-= clk.TGD1;
    sv->ts= time-sv->tsdt;
    return 0;
}

int EphemBst::ReadFromFile(string fileName, SvAll &svs) {
    return 0;
}

void *EphemBst::DownLoadThread(void *_gnss) {
    return nullptr;
}

bool EphemBst::Available(GnssTime time) {
    bool result;
    result = framOk[0]*framOk[1]*framOk[2];
    if(sv->isBeiDouGEO){
        for(int i = 3;i<10;i++){
            result*=framOk[i];
        }
    }
    if(!result){
        sprintf(sv->tip,"LackFrame");
    }
    if(SatH1){
        sv->healthy= false;
        sprintf(sv->tip,"noHeal");
    }
    return result;
}


int EphemBst::DecodeSubFrame(uint32_t *dwrds) {
    switch (sv->type){
        case SYS_GPS:
            DecodeGps(dwrds);
            break;
        case SYS_BDS:
            if(sv->isBeiDouGEO)DecodeBdsD2Frame1(dwrds);
            else DecodeBdsD1(dwrds);
            break;
    }
    return 0;
}

int EphemBst::DecodeGps(uint32_t *dwrds) {
//    printf("Decode Gps \n" );
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
    framOk[frame-1] = 1;
    printf(" Frame SYS_GPS  frame:%d",frame);

    uint32_t L2,PCodeState;
    switch(frame){
        case 1:
            clk.WN = Read1Word(dwrds[2],10,2);
            L2 = Read1Word(dwrds[2],2,12);
            URAI = Read1Word(dwrds[2],4,14);
            SatH1 = Read1Word(dwrds[2],6,18);
            clk.AODC = Read2Word(dwrds[3],2,24,dwrds[7],8,2);
            PCodeState = Read1Word(dwrds[3],1,2);
            clk.TGD1 = (int32_t)Read1Word(dwrds[6],8,18,true)*pow(2,-31);
            clk.toc = Read1Word(dwrds[7],16,10)*pow(2,4);
            clk.a2 = (int32_t)Read1Word(dwrds[8],8,2,true)*pow(2,-55);
            clk.a1 = (int32_t)Read1Word(dwrds[8],16,10,true)*pow(2,-43);
            clk.a0 = (int32_t)Read1Word(dwrds[9],22,2,true)*pow(2,-31);
            break;
        case 2:
            //todo:IODE judge?
            clk.AODE = Read1Word(dwrds[2],8,2);
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
            clk.AODE = Read1Word(dwrds[9],8,2);
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
    return 0;
}

int EphemBst::DecodeBdsD1(uint32_t *dwrds) {
//    printf("Decode BDSD1 \n" );
    if(1810!=Read1Word(dwrds[0],11,2))
        return false;
    int frame = Read1Word(dwrds[0],3,17);
    framOk[frame-1] = 1;
    clk.SOW = Read2Word(dwrds[0],8,20,dwrds[1],12,2);
//    printf(" Frame BeidouD1 svid:%d,frame:,%d\n",svId,frame);
    switch (frame){
        case 1:
            SatH1 = Read1Word(dwrds[1],1,14);
            URAI = Read1Word(dwrds[1],4,20);
            if(URAI)printf("\n\n\nUARI not ok = %d\n\n\n",URAI);
            clk.WN = Read1Word(dwrds[2],13,2);
            ion->a0 = ((int32_t) Read1Word(dwrds[4],8,8,true))*pow(2,-30);
            ion->a1 = ((int32_t) Read1Word(dwrds[4],8,16,true))*pow(2,-27)/GPS_PI;
            ion->a2 = ((int32_t) Read1Word(dwrds[5],8,2,true))*pow(2,-24)/GPS_PI2;
            ion->a3 = ((int32_t) Read1Word(dwrds[5],8,10,true))*pow(2,-24)/GPS_PI3;
            ion->b0 = ((int32_t) Read2Word(dwrds[5],6,18,dwrds[6],2,2,true))*pow(2,11);
            ion->b1 = ((int32_t) Read1Word(dwrds[6],8,4,true))*pow(2,14)/GPS_PI;
            ion->b2 = ((int32_t) Read1Word(dwrds[6],8,12,true))*pow(2,16)/GPS_PI2;
            ion->b3 = ((int32_t) Read2Word(dwrds[6],4,20,dwrds[7],4,2,true))*pow(2,16)/GPS_PI3;

            clk.AODC = Read1Word(dwrds[1],5,15);

            clk.toc = Read2Word(dwrds[2],9,15,dwrds[3],8,2)*8;
            clk.a0 = (int32_t)Read2Word(dwrds[7],7,17,dwrds[8],17,2,true)*pow(2,-33);
            clk.a1 = (int32_t)Read2Word(dwrds[8],5,19,dwrds[9],17,2,true)*pow(2,-50);
            clk.a2 = (int32_t)Read1Word(dwrds[7],11,6,true)*pow(2,-66);

            clk.TGD1 = (int32_t)Read1Word(dwrds[3],10,10,true)*1e-10;
            clk.TGD2 = (int32_t)Read2Word(dwrds[3],4,20,dwrds[4],6,2,true)*1e-10;
            clk.AODE = Read1Word(dwrds[9],5,19);
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
    return 0;
}

int EphemBst::DecodeBdsD2Frame1(uint32_t *dwrds) {
//    printf("Decode BDSD2 \n" );
    if(1810!=Read1Word(dwrds[0],11,2))return -1;
    int frame = Read1Word(dwrds[0],3,17);
    if(1!=frame)return -1;
    clk.SOW = Read2Word(dwrds[0],8,20,dwrds[1],12,2);
    int Pnum1 = Read1Word(dwrds[1],4,14);
    framOk[Pnum1-1] = 1;

//    printf(" Frame BeidouD2 svid:%d,frame1,page:,%d",svId,Pnum1);
    switch(Pnum1){
        case 1:
            SatH1 = Read1Word(dwrds[1],1,18);
            clk.AODC = Read1Word(dwrds[1],5,19);
            URAI = Read1Word(dwrds[2],4,2);
            if(URAI)printf("\n\n\nUARI not ok = %d\n\n\n",URAI);
            clk.WN = Read1Word(dwrds[2],13,6);
            clk.toc = Read2Word(dwrds[2],5,19,dwrds[3],12,2)*8;
            clk.TGD1 = (int32_t)Read1Word(dwrds[3],10,14,true)*1e-10;
            clk.TGD2 = (int32_t)Read1Word(dwrds[4],10,2,true)*1e-10;
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
            clk.a0 = (int32_t)Read2Word(dwrds[3],12,12,dwrds[4],12,2,true)*pow(2,-33);
            clk.a1High = Read1Word(dwrds[4],4,14,true)<<18;
            clk.a1 = (int32_t)(clk.a1High|clk.a1Low)*pow(2,-50);
            break;
        case 4:
            clk.a1Low = Read2Word(dwrds[1],6,18,dwrds[2],12,2);
            clk.a1 = (int32_t)(clk.a1High|clk.a1Low)*pow(2,-50);
            clk.a2 = (int32_t)Read2Word(dwrds[2],10,14,dwrds[3],1,2,true)*pow(2,-66);
            clk.AODE = Read1Word(dwrds[3],5,3);
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
    return 0;
}

uint32_t EphemBst::Read1Word(uint32_t word, int length, int head, bool isInt) {
    if(isInt)
        return uint32_t (((int32_t)word)<<head>>(32-length));
    else
        return word<<head>>(32-length);
}

uint32_t
EphemBst::Read2Word(uint32_t word0, int length0, int head0, uint32_t word1, int length1, int head1, bool isInt) {
    uint32_t high,low;
    if(isInt)
        high = uint32_t (((int32_t)word0)<<head0>>(32-length0)<<length1);
    else
        high = word0<<head0>>(32-length0)<<length1;
    low = word1<<head1>>(32-length1);
    return  high|low;
}

uint32_t
EphemBst::Read3Word(uint32_t word0, int length0, int head0, uint32_t word1, int length1, int head1, uint32_t word2,
                    int length2, int head2, bool isInt) {
    uint32_t high,low;
    high = Read2Word(word0,length0,head0,word1,length1,head1,isInt)<<length2;
    low = word2<<head2>>(32-length2);
    return  high|low;
}
