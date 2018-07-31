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
    char loginstring[512] = "GET /%s HTTP/1.0\r\n"
                            "User-Agent: NTRIP GNSSInternetRadio/1.4.10\r\n"
                            "Accept: */*\r\n"
                            "Connection: close\r\n"
                            "Authorization: Basic cXhhb3lhMDAxOmEwMjk3Yzg=\r\n\r\n"; //base64 code: username:password
    char buffer[512];
    sprintf(buffer, loginstring, rtk_protocol.c_str());

    char loginReturn[512] = {0};
    send(sock_, buffer, strlen(buffer), 0);        // 向服务器发送信息
    recv(sock_, loginReturn, sizeof(loginReturn), 0);    // 接收数据
    printf("login = %s", loginReturn);
    if (0 == memcmp("ICY 200 OK", loginReturn, 10)) {
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
            //printf("get %d bytes from cors/qx\n", (int)recvLength);
//            if(sp_->write_some(boost::asio::buffer(bufferRecv, recvLength)) != recvLength) {
//                printf("failed to write data to serial port\n");
//            }
            //todo:paraseRTK;
//            fwrite(bufferRecv, recvLength, 1, fp);
            for(int i = 0;i<recvLength;i++){
                if(0xd3==bufferRecv[i]){
                    uint32_t messageLength = NetToHost32(bufferRecv+i,14,10);
                    uint32_t checkSumGet = NetToHost32(bufferRecv+i+3+messageLength,0,24);
                    crc24Q.process_bytes(bufferRecv,messageLength+3);
                    if(checkSumGet==crc24Q.checksum()){
                        memcpy(buferRTK,bufferRecv+i+3,messageLength);
                        uint32_t type = NetToHost32(buferRTK,0,12);
                        switch (type){
                            case 1005:
                                ParaseRtk32_1005(buferRTK);
                                break;
                            case 1074:
                                break;
                            case 1084:
                                break;
                            case 1124:
                                break;
                            default:
                                break;
                        }

                        i+=messageLength;
                    }
                }
            }
            ParaseMSM4(bufferRecv);
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


int NtripRTK::ParaseMSM4(const char *bufferRTK) {
    //first Parase MSM message Header:
    //type:
    uint16_t type = *(uint16_t*)bufferRTK>>4;
    uint16_t RefID = *(uint16_t*)(bufferRTK+1)<<4>>4;
    uint32_t rtkTime = *(uint32_t*)(bufferRTK+3)>>2;
    uint8_t multiMsg = *(uint8_t*)(bufferRTK+6)<<6>>7;
    uint16_t IODS = *(uint16_t*)(bufferRTK+6)<<7>>13;
    uint8_t clockCorrect = *(uint8_t*)(bufferRTK+8)<<1>>6;
    uint8_t clockExtern = *(uint8_t*)(bufferRTK+8)<<3>>6;
    uint8_t smoothType = *(uint8_t*)(bufferRTK+8)<<5>>7;
    uint16_t smoothRange = *(uint16_t*)(bufferRTK+8)<<6>>13;
    uint64_t satMask = (*(uint64_t*)(bufferRTK+9)>>56<<57)||(*(uint64_t*)(bufferRTK+10)>>7);
    uint32_t sigMask = (*(uint32_t*)(bufferRTK+17)>>24<<25)||(*(uint32_t*)(bufferRTK+18)>>7);
}



uint32_t NtripRTK::NetToHost32(char *begin, int head, int length) {
    char data[4];
    data[0] = begin[3];
    data[1] = begin[2];
    data[2] = begin[1];
    data[3] = begin[0];
    return  *(uint32_t*)data<<head>>(32-length);
}

//todo:test 1005

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