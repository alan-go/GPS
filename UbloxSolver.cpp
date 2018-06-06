#include "UbloxSolver.h"
using namespace std;


UbloxSolver::UbloxSolver(){}
UbloxSolver::~UbloxSolver(){}

bool UbloxSolver::ParseRawData(char *message) {
    cout<<"start solving rawdata"<<endl;
    char* playload = message + 6;

    char* temp = playload;
    rcvtow = *(double*)temp;
    temp = playload + 11;
    numMeas = *(u_int8_t*)temp;
    if(0==numMeas)return false;

    vector<int> GPSIndexTemp,BeiDouIndexTemp;
    for(u_int8_t n = 0;n<numMeas;n++){
        SvInfo* svTemp;
        int n32 = n*32;
        uint8_t gnssId = *(uint8_t *)(playload+36+n32);
        uint8_t svId = *(uint8_t *)(playload+37+n32);
        switch (gnssId){
            case 0:
                GPSIndexTemp.push_back(svId);
                svTemp = &(GPSSVs[svId]);
            case 3:
                BeiDouIndexTemp.push_back(svId);
                svTemp = &(BeiDouSVs[svId]);
        }
        svTemp->prMes = *(double*)(playload+16+n32);
        svTemp->cpMes = *(double*)(playload+24+n32);
        svTemp->doMes = *(float *)(playload+32+n32);

    }
    GPSIndex = GPSIndexTemp;
    BeiDouIndex = BeiDouIndexTemp;

    positionMtx.lock();
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
    //todo:
    positionMtx.unlock();
}

bool UbloxSolver::DecodeBeiDouBroadcastD1(uint32_t *dwrds, SvInfo *sv) {
    //todo:
    uint32_t bit32,bit32x2[2];
    uint64_t bit64;
    //word0:
    if(1810!=Read1Word(dwrds[0],11))
        return false;
    int page = Read1Word(dwrds[0],3,17);
    switch (page){
        case 1:
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
            sv->orbit.toeF2 = Read1Word(dwrds[9],2,22)<<15<<3;
            sv->orbit.toe = sv->orbit.toeF2|sv->orbit.toeF3;
            sv->orbit.sq_a = Read2Word(dwrds+8,12,12,20)*pow(2,-19);
            sv->orbit.e = Read2Word(dwrds+4,10,14,22)*pow(2,-33);
            sv->orbit.dtn = (int32_t)Read2Word(dwrds+1,10,14,6,2,true)*pow(2,-43);
            sv->orbit.M0 = (int32_t)Read2Word(dwrds+3,20,4,12,2,true)*pow(2,-31);

        case 3:
            sv->orbit.toeF3 = Read2Word(dwrds+2,10,14,5)<<3;
            sv->orbit.toe = sv->orbit.toeF2|sv->orbit.toeF3;
            sv->orbit.omega = (int32_t)Read2Word(dwrds+8,11,13,21,2,true)*pow(2,-31);
            sv->orbit.Omega0 = (int32_t)Read2Word(dwrds+7,21,3,11,2,true)*pow(2,-31);
            sv->orbit.OmegaDot = (int32_t)Read2Word(dwrds+4,11,13,13,2,true)*pow(2,-43);
            sv->orbit.i0 = (int32_t)Read2Word(dwrds+2,17,7,15,2,true)*pow(2,-31);
            sv->orbit.IDOT = (int32_t)Read2Word(dwrds+6,13,11,1,2,true)*pow(2,-43);

    }

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