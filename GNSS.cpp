#include "GNSS.h"

GNSS::GNSS() {

}

GNSS::~GNSS() {

}

int GNSS::StartGNSS(const std::string &serial_port, const unsigned int baudRate) {
    //    rtk_protocol_ = rtk_protocol;
    serialDataManager.serialPort_ = serial_port;
    serialDataManager.baudRate = baudRate;
    serialDataManager.stopCapture = false;
    //todo ntrip false;

    //todo:NtripLogin
//    if(!NtripLogin(rtk_protocol_)) {
//        printf("cannot login ntrip server\n");
//        exit(0);
//    }

    pthread_create(&thread2_, nullptr, ThreadAdapterQianXun, this);
    sleep(2);
    pthread_create(&thread1_, nullptr, ThreadAdapterGNSS, &serialDataManager);
    sleep(2);
    return 1;
}

void* GNSS::ThreadAdapterGNSS(void *__this) {
    auto _this= ( SerialData* ) __this;
    _this->StartCapture(_this->serialPort_,_this->baudRate);
    return nullptr;
}

void* GNSS::ThreadAdapterQianXun(void *__this) {
    auto _this= ( GNSS* ) __this;

//    todo:
        // _this->RecvThread();
    return nullptr;
}