#include "NtripRTK.h"
#include "GNSS.h"

NtripRTK::NtripRTK():crc24Q( 0x864cfb, 0, 0, false, false ) ,refStationId(0),isPhysicalStation(0) ,
singleReceiver(0),ITRFyear(0),supportBeiDou(0),supportGalileo(0),supportGLONASS(0),supportGPS(0),quaterCycle(0){
    ECEF_XYZ<<0,0,0;
}

NtripRTK::~NtripRTK() {

}

bool NtripRTK::NtripLogin(const std::string &rtk_protocol) {
    sock_ = socket(AF_INET, SOCK_STREAM, 0);
    if (sock_ < 0) {
        perror("socket");
        return false;
    }

    // 设置服务器地址结构体
    struct sockaddr_in server_addr;
    bzero(&server_addr, sizeof(server_addr)); // 初始化服务器地址
    server_addr.sin_family = AF_INET;   // IPv4
    server_addr.sin_port = htons(port_); // 端口
    inet_pton(AF_INET, serverIP_.c_str(), &server_addr.sin_addr);   // ip

    // 连接服务器
    int err_log = connect(sock_, (struct sockaddr *) &server_addr, sizeof(server_addr));
    if (err_log != 0) {
        perror("connect");
        close(sock_);
        return false;
    }

    //登录请求
    char loginstring[512] = "GET /RTCM32_GGB HTTP/1.0\r\n"
                            "User-Agent: NTRIP GNSSInternetRadio/1.4.10\r\n"
                            "Accept: */*\r\n"
                            "Connection: close\r\n"
//                            "Authorization: Basic cXhhb3lhMDAxOmEwMjk3Yzg=\r\n\r\n"; //base64 code: username:password
                            "Authorization: Basic cXhhb3lhMDAxOmEwMjk3Yzg=\r\n\r\n"; //findMM
    char buffer[512];
    sprintf(buffer, loginstring, rtk_protocol.c_str());

    char loginReturn[512] = {0};
    send(sock_, buffer, strlen(buffer), 0);        // 向服务器发送信息
    recv(sock_, loginReturn, sizeof(loginReturn), 0);    // 接收数据
    printf("loginReturn = %s", loginReturn);
    if (0 == memcmp("ICY 200 OK", loginReturn, 10)) {
        printf("Login success.\n");
        return true;
    }
    else {
        printf("\nntrip login failed.\n");
        close(sock_);
        return false;
    }
    return true;
}


void NtripRTK::RecvThread() {
    char bufferRecv[1024];
    char buferRTK[256];
    printf("in recv thread\n");
    FILE *fp = fopen("cmr.bin", "wb");
    while(!stopRTK) {
        auto recvLength = recv(sock_, bufferRecv, sizeof(bufferRecv), 0);
        if(recvLength < 0) {
            close(sock_);
            printf("socket recv error! try to login again");
            gnss->StopGNSS();
            NtripLogin(gnss->rtk_protocol_);
        }
        if(recvLength >= 0) {
            for(int i = 0;i<recvLength;i++){
                if(0xd3==(u_char)bufferRecv[i]&&0==(bufferRecv[i+1]>>2)){
                    printf("get 0xd3,i = %d.\n",i);
                    uint32_t messageLength = NetToHost32(bufferRecv+i,14,10);
                    printf("messageLength = %d\n",messageLength);
                    uint32_t checkSumGet = NetToHost32(bufferRecv+i+3+messageLength,0,24);
                    crc24Q.reset();
                    crc24Q.process_bytes(bufferRecv+i,messageLength+3);
                    if(checkSumGet==crc24Q.checksum()){
                        printf("CheckSum OK.\n");
                        memcpy(buferRTK,bufferRecv+i+3,messageLength);
                        uint32_t type = NetToHost32(buferRTK,0,12);
                        switch (type){
                            case 1005:
                                ParaseRtk32_1005(buferRTK);
                                break;
                            case 1074:
                                ParaseMSM4(buferRTK, SV::GPS);
                                break;
                            case 1084:
                                break;
                            case 1124:
                                ParaseMSM4(buferRTK, SV::BeiDou);
                                break;
                            default:
                                break;
                        }

                        i+=messageLength;
                    } else{
                        printf("CheckSum failed.%d,%d\n",checkSumGet,crc24Q.checksum());
                    }
                }
            }
        }
    }
    fclose(fp);
}

int NtripRTK::SentGGA(const char *bufferGGA, int length) {
    //todo UpdateGGA;
    int ret = send(sock_, bufferGGA, length, 0);
    if(ret != length) {
        printf("socket send error! try to login again\n");
        gnss->StopGNSS();
        NtripLogin(gnss->rtk_protocol_);
    }
}

void NtripRTK::UpdateGGA() {}


int NtripRTK::ParaseMSM4(char *bufferRTK, SV::SvType type) {
    //first Parase MSM message Header:
    uint32_t RefID = NetToHost32(bufferRTK,12,12);
    uint32_t rtkTime = NetToHost32(bufferRTK+3,0,30);
    uint8_t multiMsg = *(uint8_t*)(bufferRTK+6)<<6>>7;
    uint32_t IODS = NetToHost32(bufferRTK+6,7,3);
    uint8_t clockCorrect = *(uint8_t*)(bufferRTK+8)<<1>>6;
    uint8_t clockExtern = *(uint8_t*)(bufferRTK+8)<<3>>6;
    uint8_t smoothType = *(uint8_t*)(bufferRTK+8)<<5>>7;
    uint32_t smoothRange = NetToHost32(bufferRTK+8,6,3);

    vector<MSM4data*> satsData;
    vector<MSM4Cell*> cells;
//    vector<SV::SignalData*>sigData;
    int nSat = 0, nSig = 0, nCell = 0;

    //sats里面是卫星编号，对应手册
    vector<int> sats,sigs;


    //sat
    char *b = bufferRTK+9;
    int i = 1;
    u_char temp = (u_char)b[0];
    //卫星掩膜
    for(int id=0;id<64;id++,i++){
        if(8==i){
            i=0;
            b++;
            temp = (u_char)b[0];
        }
        if(temp&(128>>i)){
            sats.push_back(id);
            nSat++;
            printf("Sat:%d\n",id+1);
            MSM4data *tempData = new MSM4data;
            tempData->rtktime = rtkTime;
            satsData.push_back(tempData);
            AddRtkRecord(tempData,type,id);
        }
    }
    //signal信号掩膜
    for(int id=0;id<32;id++,i++){
        if(8==i){
            i=0;
            b++;
            temp = (u_char)b[0];
        }
        if(temp&(128>>i)){
            sigs.push_back(id+1);
            nSig++;
            printf("Sig:%d\n",id+1);
        }
    }
    //cell
    for(int isig = 0;isig<nSig;isig++){
        for(int isat = 0;isat<nSat;isat++){
            if(8==i){
                i=0;
                b++;
                temp = (u_char)b[0];
            }
            if(temp&(128>>i)){
                //todo
                nCell++;
//                sigData.push_back(svHead[sats[isat]].SignalTable(sigs[isig]));
                int m = sats[isat];
                int n = sigs[isig];
                printf("m,n=%d,%d\n",m,n);
                MSM4data *dataTemp = GetRtkRecord(m,0,type);
                cells.push_back(dataTemp->sigData+n);
                dataTemp->sigs.push_back(n);
//                SV * temp = gnss->svsManager.SatTable(type,m);
//                sigData.push_back(temp->SignalTable(n));
            }
            i++;
        }
    }

    //parase Sat data
    for(int n = 0;n<nSat;n++){
        satsData[n]->df397 = double(NetToHost32(b,i,8));
//        *(df397[n]) = double(NetToHost32(b,i,8));
        printf("df397=%lf\n",satsData[n]->df397);
        b++;
    }
    for(int n = 0;n<nSat;n++){
        if(i>=8){
            i-=8;
            b++;
        }
        satsData[n]->df398 = double(NetToHost32(b,i,10));

//        *(df398[n]) = double(NetToHost32(b,i,10));
        printf("df398=%lf\n",satsData[n]->df398);
        b++;
        i+=2;
    }
    //parase sig data
    for(int n = 0; n< nCell; n++){
        if(i>=8){
            i-=8;
            b++;
        }
        cells[n]->df400 = double(NetToHost32(b,i,15));
//        sigData[n]->df400 = double(NetToHost32(b,i,15));
        printf("df400=%lf\n",cells[n]->df400);
        b++;
        i+=7;
    }
    for(int n = 0; n< nCell; n++){
        if(i>=8){
            i-=8;
            b++;
        }
        cells[n]->df401 = double(NetToHost32(b,i,22));
        printf("df401=%lf\n",cells[n]->df401);

        b+=2;
        i+=6;
    }
    for(int n = 0; n< nCell; n++){
        if(i>=8){
            i-=8;
            b++;
        }
        cells[n]->df402 = NetToHost32(b,i,4);
        printf("df402=%d\n",cells[n]->df402);

        i+=4;
    }
    for(int n = 0; n< nCell; n++){
        if(i>=8){
            i-=8;
            b++;
        }
        cells[n]->df420 = NetToHost32(b,i,1);
        i++;
    }
    for(int n = 0; n< nCell; n++){
        if(i>=8){
            i-=8;
            b++;
        }
        cells[n]->df403 = NetToHost32(b,i,6);
        printf("df403=%d\n",cells[n]->df403);

        i+=6;
    }

//    GNSS * g2 = gnss;
//    printf("df400=%lf\n",g2->svsManager.svGpss[11].L1_1C.df400);
//    printf("df400=%lf\n",gnss->svsManager.svGpss[12].L1_1C.df400);
//    printf("df400=%lf\n",g2->svsManager.svGpss[13].L1_1C.df400);
//    printf("df400=%lf\n",g2->svsManager.svGpss[14].L1_1C.df400);
//    printf("df397=%lf\n",g2->svsManager.svGpss[11].df397);
//    printf("df397=%lf\n",g2->svsManager.svGpss[13].df397);
//    for(int i=0;i<nSat;i++) {
//        int ind = sats[i];
//        printf("ind = %d\n",ind);
//        SV *temp1 = gnss->svsManager.SatTable(SV::GPS, ind);
//        double pp = temp1->df397 + temp1->df398 / 1024;
//        double c1 = pp + (temp1->SignalTable(2)->df400) * pow(2,-24);
//        double l1 = pp + temp1->SignalTable(2)->df401 / 1024 / 1024 / 512;
//        printf("c1=%lf,l1=%lf\n", pp * Light_speed, l1 * Light_speed);
//    }
}



uint32_t NtripRTK::NetToHost32(char *begin, int head, int length) {
    char data[4];
    data[0] = begin[3];
    data[1] = begin[2];
    data[2] = begin[1];
    data[3] = begin[0];
    return  *(uint32_t*)data<<head>>(32-length);
}

int NtripRTK::ParaseRtk32_1005(char *buffer) {
    refStationId = NetToHost32(buffer,12,12);
    ITRFyear = NetToHost32(buffer+3,0,6);
    if(1==NetToHost32(buffer+3,6,1))supportGPS = true;
    if(1==NetToHost32(buffer+3,7,1))supportGLONASS = true;
    if(1==NetToHost32(buffer+4,0,1))supportGalileo = true;
    if(0==NetToHost32(buffer+4,1,1))isPhysicalStation = true;

    int64_t ARP_X = *(int64_t*)(buffer+4)<<58>>26;
    uint32_t ARP_Xl = NetToHost32(buffer+5,0,32);
    ARP_X = ARP_X | ARP_Xl;

    if(1==NetToHost32(buffer+9,1,1))supportBeiDou = true;
    int64_t ARP_Y = *(int64_t*)(buffer+9)<<58>>26;
    uint32_t ARP_Yl = NetToHost32(buffer+10,0,32);
    ARP_Y = ARP_Y | ARP_Yl;

    quaterCycle = NetToHost32(buffer+14,0,2);
    int64_t ARP_Z = *(int64_t*)(buffer+14)<<58>>26;
    uint32_t ARP_Zl = NetToHost32(buffer+15,0,32);
    ARP_Z = ARP_Z | ARP_Zl;

    ECEF_XYZ<<ARP_X,ARP_Y,ARP_Z;
    ECEF_XYZ*=1e-4;
}

int NtripRTK::TestParase(char *bufferRecv,int recvLength) {
    char buferRTK[256];

    for(int i = 0;i<recvLength;i++){
        if(0xd3==(u_char)bufferRecv[i]&&0==(bufferRecv[i+1]>>2)){
            printf("get 0xd3,i = %d.\n",i);
            uint32_t messageLength = NetToHost32(bufferRecv+i,14,10);
            printf("messageLength = %d\n",messageLength);
            uint32_t checkSumGet = NetToHost32(bufferRecv+i+3+messageLength,0,24);
            crc24Q.reset();
            crc24Q.process_bytes(bufferRecv+i,messageLength+3);
            if(checkSumGet==crc24Q.checksum()){
                printf("CheckSum OK.\n");
                memcpy(buferRTK,bufferRecv+i+3,messageLength);
                uint32_t type = NetToHost32(buferRTK,0,12);
                switch (type){
                    case 1005:
                        ParaseRtk32_1005(buferRTK);
                        break;
                    case 1074:
                        ParaseMSM4(buferRTK, SV::GPS);
                        break;
                    case 1084:
                        break;
                    case 1124:
                        ParaseMSM4(buferRTK, SV::BeiDou);
                        break;
                    default:
                        break;
                }

                i+=messageLength;
            } else{
                printf("CheckSum failed.%d,%d\n",checkSumGet,crc24Q.checksum());
            }
        }
    }
//    ParaseMSM4(bufferRecv);
}

int NtripRTK::AddRtkRecord(MSM4data *data, SV::SvType type, int id) {
    //maximum:5
    int maxNumber = 5;
    vector<MSM4data*> *temp;
    switch (type){
        case SV::GPS:
            temp = rtkDataGps+id;
//            rtkDataGps[id].insert(rtkDataGps[id].begin(),data);
//            if(rtkDataGps[id].size()>maxNumber)rtkDataGps[id].pop_back();
            break;
        case SV::BeiDou:
            temp = rtkDataBeiDou+id;
//            rtkDataBeiDou[id].insert(rtkDataBeiDou[id].begin(),data);
//            if(rtkDataBeiDou[id].size()>maxNumber)rtkDataBeiDou[id].pop_back();
            break;
        default:
            break;
    }
    temp->insert(temp->begin(),data);
    if(temp->size()>maxNumber){
        MSM4data* endp = *temp->end();
        temp->pop_back();
        delete(endp);
    }

    return 0;
}

MSM4data* NtripRTK::GetRtkRecord(int satInd, int timeInd, SV::SvType type) {
    switch (type){
        case SV::GPS:
            return rtkDataGps[satInd][timeInd];
        case SV::BeiDou:
            return rtkDataBeiDou[satInd][timeInd];
        default:
            return nullptr;

    }
}

