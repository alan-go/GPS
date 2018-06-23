#include "UbloxSolver.h"

using namespace std;

SvInfo::SvInfo(){
    page1OK = page2OK = page3OK = false;
}
SvInfo::~SvInfo(){}

bool SvInfo::CalcuTime(double rcvtow) {
    ts = rcvtow - prMes/M_PI;
    double t = ts-toc;
    //此时的钟差是没有考虑相对论效应和 TGD的
    for(int i = 0;i<5;i++){
        tsDelta = a0+a1*t+a2*t*t;
        t -=tsDelta;
    }
    tsReal = ts - tsDelta;
}

 bool SvInfo::CalcuECEF(double rcvtow) {
    CalcuTime(rcvtow);
    double A = orbit.sq_a*orbit.sq_a;
    double n0 = sq_M_miu/(orbit.sq_a*orbit.sq_a*orbit.sq_a);
    double tk = tsReal - orbit.toe;
    double n = n0 + orbit.dtn;
    double Mk = orbit.M0 + n*tk;
    double Ek = Mk,EkOld = Ek+1;
    while(Ek-EkOld>1e-8){
        EkOld = Ek;
        Ek = EkOld-(EkOld-orbit.e*sin(EkOld)-Mk)/(1-orbit.e*cos(EkOld));
    }
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
    double Omegak = orbit.Omega0 + (orbit.OmegaDot-Omega_e)*tk - Omega_e*orbit.toe;
    position(0) = xk*cos(Omegak) - yk*cos(ik)*sin(Omegak);
    position(1) = xk*sin(Omegak) + yk*cos(ik)*cos(Omegak);
    position(2) = yk*sin(ik);
}

UbloxSolver::UbloxSolver(){
    rxyz = Eigen::Vector4d::Zero();
}

UbloxSolver::~UbloxSolver(){}


bool UbloxSolver::ParseRawData(char *message) {
    cout<<"start solving rawdata"<<endl;
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
                svTemp = &(GPSSVs[svId]);
                svTemp->type = SvInfo::GPS;
            case 3:
                svTemp = &(BeiDouSVs[svId]);
                svTemp->type = SvInfo::BeiDou;
        }
        svTemp->svId = svId;
        visibleSvs.push_back(svTemp);
        svTemp->prMes = *(double*)(playload+16+n32);
        svTemp->cpMes = *(double*)(playload+24+n32);
        svTemp->doMes = *(float *)(playload+32+n32);

    }

    //将星历拷贝一份防止计算时被改.
    vector<SvInfo>().swap(SvsForCalcu);
    for(int i=0;i<numMeas;i++)
    {
        SvInfo sv = *visibleSvs[i];
        if(sv.page123OK)SvsForCalcu.push_back(sv);
    }
    if(SvsForCalcu.size()<5)
        return false;

    thread PositionThread(LaunchPositionThread,this);
    PositionThread.detach();
    return true;
}

bool UbloxSolver::ParseBstSubFrame(char *message) {
    cout<<"Update subframe"<<endl;
    char* playload = message + 6;

    uint8_t gnssId = *(uint8_t*)(playload);
    uint8_t svId = *(uint8_t*)(playload + 1);
    uint8_t numWords = *(uint8_t*)(playload+4);

    char* tmp = playload+8;

    if(10==numWords){
        uint32_t dwrds[10];
        for(int i=0;i<10;i++)   dwrds[i] = *(uint32_t*)(tmp+4*i);
        switch (gnssId){
            case 0:
                DecodeGpsBroadcast(dwrds,&GPSSVs[svId]);
            case 3:
                if(svId>5)
                    DecodeBeiDouBroadcastD1(dwrds,&BeiDouSVs[svId]);
                else
                    DecodeBeiDouBroadcastD2(dwrds,&BeiDouSVs[svId]);
        }
    }

}

bool UbloxSolver::solvePosition() {
    if(isCalculating)
        return false;
    //todo:
    //选星？
    int N = SvsForCalcu.size();
    rxyzOld = rxyz;
    Eigen::Vector4d dtxyzt = Eigen::Vector4d::Ones();
    Eigen::MatrixXd pc(N,1);
    Eigen::MatrixXd b(N,1);

    for(int i = 0;i< N;i++){
        SvInfo* sv= &SvsForCalcu[i];
        sv->CalcuECEF(rcvtow);
        //todo: calcu I,T,dtu
        double pci = sv->prMes + sv->dts - sv->I - sv->T;
        pc(i) = pci;
    }

    int numCalcu = 0;
    while (dtxyzt.squaredNorm()>0.1){
        numCalcu++;
        Eigen::MatrixXd G(N,4);
        for(int i = 0;i< N;i++){
            SvInfo* sv= &SvsForCalcu[i];
            double r = (sv->position-rxyz.head(3)).squaredNorm();
            b(i) = pc(i)-r-rxyz(3);
            G(i,0) = (rxyz(0)-sv->position(0))/r;
            G(i,1) = (rxyz(1)-sv->position(1))/r;
            G(i,2) = (rxyz(2)-sv->position(2))/r;
            G(i,3) = 1;
        }
        Eigen::Matrix GT = G.transpose();
        dtxyzt = ((GT*G).inverse())*GT*b;
        if(numCalcu>30)return false;
    }

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
    switch (page){
        case 1:
            sv->page1OK = true;
            sv->SOW = Read1Word(dwrds[0],8,20);
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
        case 3:
            sv->page3OK = true;
            sv->orbit.toeF3 = Read2Word(dwrds+2,10,14,5)<<3;
            sv->orbit.toe = sv->orbit.toeF2|sv->orbit.toeF3;
            sv->orbit.omega = (int32_t)Read2Word(dwrds+8,11,13,21,2,true)*pow(2,-31);
            sv->orbit.Omega0 = (int32_t)Read2Word(dwrds+7,21,3,11,2,true)*pow(2,-31);
            sv->orbit.OmegaDot = (int32_t)Read2Word(dwrds+4,11,13,13,2,true)*pow(2,-43);
            sv->orbit.i0 = (int32_t)Read2Word(dwrds+2,17,7,15,2,true)*pow(2,-31);
            sv->orbit.IDOT = (int32_t)Read2Word(dwrds+6,13,11,1,2,true)*pow(2,-43);
            sv->orbit.Cic = (int32_t)Read2Word(dwrds+3,7,17,11,2,true)*pow(2,-31);
            sv->orbit.Cis = (int32_t)Read2Word(dwrds+5,9,15,9,2,true)*pow(2,-31);
    }
    sv->page123OK = sv->page1OK && sv->page3OK && sv->page3OK;
}

bool UbloxSolver::DecodeBeiDouBroadcastD2(uint32_t *dwrds, SvInfo *sv) {

}

bool UbloxSolver::DecodeGpsBroadcast(uint32_t *dwrds, SvInfo *sv) {
    //todo:
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