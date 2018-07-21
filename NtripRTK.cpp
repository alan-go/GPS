#include "NtripRTK.h"
#include "GNSS.h"

NtripRTK::NtripRTK() {

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
    char bufferRTK[256];
    printf("in recv thread\n");
    FILE *fp = fopen("cmr.bin", "wb");
    while(!stopRTK) {
        auto rtkLength = recv(sock_, bufferRTK, sizeof(bufferRTK), 0);
        if(rtkLength < 0) {
            close(sock_);
            printf("socket recv error! try to login again");
            gnss->StopGNSS();
            NtripLogin(gnss->rtk_protocol_);
        }
        if(rtkLength >= 0) {
            //printf("get %d bytes from cors/qx\n", (int)rtkLength);
//            if(sp_->write_some(boost::asio::buffer(bufferRTK, rtkLength)) != rtkLength) {
//                printf("failed to write data to serial port\n");
//            }
            //todo:paraseRTK;
            fwrite(bufferRTK, rtkLength, 1, fp);
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