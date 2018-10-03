#include "NtripRTK.h"
#include "GNSS.h"

NtripRTK::NtripRTK():crc24Q( 0x864cfb, 0, 0, false, false ) ,refStationId(0),isPhysicalStation(0), stopRTK(0) ,
singleReceiver(0),ITRFyear(0),supportBeiDou(0),supportGalileo(0),supportGLONASS(0),supportGPS(0),quaterCycle(0){
    ECEF_XYZ<<0,0,0;
    sprintf(ggaDefault,"$GPGGA,100236150.00,3958.7565643,N,11620.3448181,E,1,06,1.0,52.911,M,0.000,M,0.0,*79\r\n");
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
    printf("in recv thread\n");
    FILE *fp;
    if(gnss->logOpen)fp = fopen(saveName, "wb");
    while(!stopRTK) {
        auto recvLength = recv(sock_, bufferRecv, sizeof(bufferRecv), 0);
        if(recvLength < 0) {
            close(sock_);
            printf("socket recv error! try to login again");
            gnss->StopGNSS();
            NtripLogin(gnss->rtk_protocol_);
        }
        if(recvLength >= 0) {
            ParaseRTK(bufferRecv,recvLength);
            if(gnss->logOpen)fwrite(bufferRecv, recvLength, 1, fp);
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


int NtripRTK::ParaseMSM4(char *bufferRTK, SysType type) {
    //first Parase MSM message Header:
    uint32_t refID = NetToHost32(bufferRTK,12,12);//DF003
    double rtkTime = double(NetToHost32(bufferRTK+3,0,30))*1e-3;//bdsTime = GPSTime -14s
    uint8_t multiMsg = NetToHost32(bufferRTK+6,6,1);//DF393
    uint32_t IODS = NetToHost32(bufferRTK+6,7,3);//DF409
    uint8_t clockCorrect = NetToHost32(bufferRTK+8,1,2);//DF411
    printf("DF411 clock correct= %d\n", clockCorrect);
    uint8_t clockExtern = NetToHost32(bufferRTK+8,3,2);//DF412
    printf("DF412 clock correct= %d\n", clockExtern);
    uint8_t smoothType = NetToHost32(bufferRTK+8,5,1);//DF417
    uint32_t smoothRange = NetToHost32(bufferRTK+8,6,3);//DF418

    printf("rtktime = %f, df003=%d,df393=%d,df409=%d,df411=%d,412=%d,417=%d,df418=%d\n",
            rtkTime,refID,multiMsg,IODS,clockCorrect,clockExtern,smoothType,smoothRange);

    vector<MSM4data*> satsData;

    typedef pair<MSM4data*,MSM4Cell*> pairCell;
    vector<pairCell> cells;

//    vector<SV::SignalData*>sigData;
    int nSat = 0, nSig = 0, nCell = 0;

    //sats里面是卫星编号，对应手册
    vector<int> sats,sigs;


    //sat
    char *byteInd = bufferRTK+9;
    int headi = 1;
    u_char temp = (u_char)byteInd[0];
    //卫星掩膜
    for(int id=0;id<64;id++,headi++){
        if(8==headi){
            headi=0;
            byteInd++;
            temp = (u_char)byteInd[0];
        }
        if(temp&(128>>headi)){
            sats.push_back(id);
            nSat++;
//            printf("time = %f,,",rtkTime);
            printf("Sat:%d, ",id+1);
            MSM4data *tempData = new MSM4data;
            tempData->rtktime = rtkTime;
            tempData->refID = refID;
            satsData.push_back(tempData);
            AddRtkRecord(tempData,type,id);
        }
    }
    //signal信号掩膜
    for(int id=0;id<32;id++,headi++){
        if(8==headi){
            headi=0;
            byteInd++;
            temp = (u_char)byteInd[0];
        }
        if(temp&(128>>headi)){
            sigs.push_back(id);
            nSig++;
            printf("  [Sig:%d] ",id+1);
        }
    }
    //cell
//    printf("cell\n");
    for(int isat = 0;isat<nSat;isat++){
        for(int isig = 0;isig<nSig;isig++){
            if(8==headi){
                headi=0;
                byteInd++;
                temp = (u_char)byteInd[0];
            }
            if(temp&(128>>headi)){
                //todo
                nCell++;
                int m = sats[isat];
                int n = sigs[isig];
//                printf("m,n=%d,%d\n",m,n);
//                MSM4data *dataTemp = GetRtkRecord(m,0,type);
                MSM4data *dataTemp = satsData[isat];
                cells.push_back(pairCell(dataTemp,&(dataTemp->sigData[n])));
                dataTemp->sigs.push_back(n);
            }
            headi++;
        }
    }

    //parase Sat data
//    printf("Sat data\n");
    for(int n = 0;n<nSat;n++){
        satsData[n]->df397 = double(NetToHost32(byteInd,headi,8));
        byteInd++;
    }
    for(int n = 0;n<nSat;n++){
        if(headi>=8){
            headi-=8;
            byteInd++;
        }
        satsData[n]->df398 = double(NetToHost32(byteInd,headi,10))*pow(2,-10);
        satsData[n]->prRough = satsData[n]->df397 + satsData[n]->df398;
//        printf("prRough = %.2f\n",satsData[n]->prRough);
        byteInd++;
        headi+=2;
    }
    //parase sig data
//    printf("Sig data\n");
    for(int n = 0; n< nCell; n++){
        if(headi>=8){
            headi-=8;
            byteInd++;
        }
        cells[n].second->df400= double(int32_t(NetToHost32(byteInd,headi,15,1)))*pow(2,-24);
        cells[n].second->prMes = (cells[n].first->prRough + cells[n].second->df400)*Light_speed*1e-3;
        byteInd++;
        headi+=7;
    }
    for(int n = 0; n< nCell; n++){
        if(headi>=8){
            headi-=8;
            byteInd++;
        }
        cells[n].second->df401 = double(int32_t(NetToHost32(byteInd,headi,22,1)))*pow(2,-29);
        cells[n].second->cpMes = (cells[n].first->prRough + cells[n].second->df401)*Light_speed*1e-3;
        byteInd+=2;
        headi+=6;
    }
    for(int n = 0; n< nCell; n++){
        if(headi>=8){
            headi-=8;
            byteInd++;
        }
        cells[n].second->df402 = NetToHost32(byteInd,headi,4);
        headi+=4;
    }
    for(int n = 0; n< nCell; n++){
        if(headi>=8){
            headi-=8;
            byteInd++;
        }
        cells[n].second->df420 = NetToHost32(byteInd,headi,1);
        headi++;
    }
    for(int n = 0; n< nCell; n++){
        if(headi>=8){
            headi-=8;
            byteInd++;
        }
        cells[n].second->df403 = NetToHost32(byteInd,headi,6);
        headi+=6;
    }

//    for(int n = 0; n< nSat; n++){
//        int m = sats[n];
////                printf("m,n=%d,%d\n",m,n);
//        MSM4data *dataTemp = GetRtkRecord(m,0,type);
//        printf("m=%d,rough p = %f\n",m,dataTemp->prRough);
//    }
//    printf("1 rough = %f\n",gnss->svsManager.svBeiDous[7].rtkData[0]->prRough);
}



uint32_t NtripRTK::NetToHost32(char *begin, int head, int length,bool isInt) {
    char data[4];
    data[0] = begin[3];
    data[1] = begin[2];
    data[2] = begin[1];
    data[3] = begin[0];
    uint32_t dataUint = *(uint32_t *) data;
    int32_t dataInt = *(uint32_t *) data;
    if (isInt) {
        return uint32_t(dataInt << head >> (32 - length));
    } else{
        return  dataUint<<head>>(32-length);
    }
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
    printf("refstation Id %d\n", refStationId);
    cout<<ECEF_XYZ<<endl;
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
                        ParaseMSM4(buferRTK, SYS_GPS);
                        break;
                    case 1084:
                        break;
                    case 1124:
                        ParaseMSM4(buferRTK, SYS_BDS);
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

int NtripRTK::AddRtkRecord(MSM4data *data, SysType sys, int id) {
    //maximum:5
//    printf("ADD record\n");
    int maxNumber = 102400;
    deque<MSM4data*> *temp;
    switch (sys){
        case SYS_GPS:
            temp = &(gnss->svsManager->svGpss[id].rtkData);
            break;
        case SYS_BDS:
            temp = &(gnss->svsManager->svBeiDous[id].rtkData);
            break;
        default:
            break;
    }
    temp->push_front(data);
    if(temp->size()>maxNumber){
        MSM4data* endp = temp->back();
        temp->pop_back();
        delete(endp);
    }

//    if (temp->size() == maxNumber) {
//        for (int i = 0; i < temp->size(); ++i) {
//            double pr = (*temp)[i]->sigData[1].prMes;
//            double cp = (*temp)[i]->sigData[1].cpMes;
//            printf("time = %f, pr =  %f, cp = %f\n", (*temp)[i]->rtktime, pr, cp);
//        }
//    }

    return 0;

}

MSM4data* NtripRTK::GetRtkRecord(int satInd, int timeInd, SysType type) {
    switch (type){
        case SYS_GPS:
            return gnss->svsManager->svGpss[satInd].rtkData[timeInd];
        case SYS_BDS:
            return gnss->svsManager->svBeiDous[satInd].rtkData[timeInd];
        default:
            return nullptr;
    }
}

int NtripRTK::ParaseRTK(char *buffer, int length) {
    char buferRTK[256];
    for(int i = 0;i<length-1;i++){
        if(0xd3==(u_char)buffer[i]&&0==(buffer[i+1]>>2)){
            printf("get 0xd3,i = %d.\n",i);
            uint32_t messageLength = NetToHost32(buffer+i,14,10);
            uint32_t checkSumGet = NetToHost32(buffer+i+3+messageLength,0,24);
            crc24Q.reset();
            crc24Q.process_bytes(buffer+i,messageLength+3);
            if(checkSumGet==crc24Q.checksum()){
                printf("CheckSum OK.\n");
                memcpy(buferRTK,buffer+i+3,messageLength);
                uint32_t type = NetToHost32(buferRTK,0,12);
                printf("messageLength = %d, type = %d\n",messageLength,type);
//                        mtxData.lock();
                switch (type){
                    case 1005:
                        ParaseRtk32_1005(buferRTK);
                        break;
                    case 1074:
                        ParaseMSM4(buferRTK, SYS_GPS);
                        break;
                    case 1084:
                        break;
                    case 1124:
                        ParaseMSM4(buferRTK, SYS_BDS);
                        break;
                    default:
                        break;
                }
//                        mtxData.unlock();
                i+=messageLength;
            } else{
                printf("CheckSum failed.%d,%d\n",checkSumGet,crc24Q.checksum());
            }
        }
    }
}

