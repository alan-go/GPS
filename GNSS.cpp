#include "GNSS.h"

GNSS::GNSS() :tu(0),tuBeiDou(0),tuGps(0),useGPS(1),useBeiDou(1){
    xyz<<0,0,0;
    svsManager.gnss = this;
    serialDataManager.gnss = this;
    rtkManager.gnss = this;
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

int GNSS::StartGNSS(const std::string &fileName) {

}

int GNSS::ParseRawData(char *message) {
//    PosSolver solver(svsManager, &rtkManager, this);
    PosSolver *solver1 = new PosSolver(svsManager, &rtkManager, this);
    memcpy(solver1->raw, message, 1024);
//    solver.CalcuPosition();
    pthread_create(&threadPos, nullptr, PositionThread, solver1);
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

void* GNSS::PositionThread(void *__pos) {
    auto _pos = ( PosSolver* ) __pos;
    _pos->CalcuPosition();
}