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
    if(1810!=Read1Word(dwrds[0],2,11))
        return false;
    int page = (dwrds[0]<<17)>>29;
    switch (page){
        case 1:
            //AODC
            int AODC = (dwrds[1]<<15)>>27;

            //toc
            bit32x2[0] = (dwrds[2]<<15)>>23;
            bit32x2[1] = (dwrds[3]>>22)<<24;
            bit64 = *(uint64_t*)bit32x2;
            sv->toc = bit64>>24;
            //a0
            bit32x2[0] = (dwrds[7]<<17)>>25;
            bit32x2[1] = (dwrds[8]>>13)<<15;
            bit64 = *(uint64_t*)bit32x2;
            sv->a0 = bit64>>15;
            //a1
            bit32x2[0] = (dwrds[8]<<19)>>27;
            bit32x2[1] = (dwrds[9]>>13)<<15;
            bit64 = *(uint64_t*)bit32x2;
            sv->a1 = bit64>>15;
            //a2
            sv->a2 =  (dwrds[7]<<6)>>21;


    }

}

bool UbloxSolver::DecodeBeiDouBroadcastD2(uint32_t *dwrds, SvInfo *sv) {

}

bool UbloxSolver::DecodeGpsBroadcast(uint32_t *dwrds, SvInfo *sv) {
    //todo:
}

uint32_t UbloxSolver::Read1Word(uint32_t word, int head, int length) {
    uint32_t result = (word<<head)>>(32-length);
    return result;
}

uint64_t UbloxSolver::Read2Word(uint32_t *word, int *head, int *length) {
    uint32_t bit32x2[2];
    uint64_t result;
    bit32x2[0] = (word[0]<<head[0])>>(32-length[0]);
    bit32x2[1] = (word[1]>>(32-head[1]-length[1]))<<(32-length[1]);
    result = *(uint64_t*)bit32x2;
    result = result>>(32-length[1]);
    return result;
}

uint64_t UbloxSolver::Read3Word(uint32_t *word, int *head, int *length) {

}