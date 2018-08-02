#include "GNSS.h"

GNSS::GNSS() :tu(0),tuBeiDou(0),tuGps(0),useGPS(1),useBeiDou(1),isPositioned(false){
    xyz<<0,0,0;
    svsManager.gnss = this;
    serialDataManager.gnss = this;
    rtkManager.gnss = this;

    rtkManager.serverIP_ = "60.205.8.49";
    rtkManager.port_ = 8002;
}

GNSS::~GNSS() {

}

int GNSS::StartGNSS(const std::string &serial_port, const unsigned int baudRate) {
    //    rtk_protocol_ = rtk_protocol;
    serialDataManager.serialPort_ = serial_port;
    serialDataManager.baudRate = baudRate;
    serialDataManager.stopCapture = false;
    rtkManager.stopRTK = false;

    if(!rtkManager.NtripLogin(rtk_protocol_)) {
        printf("cannot login ntrip server\n");
        exit(0);
    }

    pthread_create(&thread2_, nullptr, ThreadAdapterQianXun, &rtkManager);
    sleep(2);
    pthread_create(&thread1_, nullptr, ThreadAdapterGNSS, &serialDataManager);
    sleep(2);
    return 1;
}

int GNSS::StartGNSS(const std::string &fileName) {

}

int GNSS::StopGNSS() {
    try {
        if (rtkManager.sock_ > 0) {
            close(rtkManager.sock_);
        }
        rtkManager.sock_ = 0;
        rtkManager.stopRTK = true;
        serialDataManager.stopCapture = true;
        if (thread2_ > 0) {//(void*)0
            pthread_join(thread2_, nullptr);
        }
        thread2_ = 0;
        if (thread1_ > 0) {//(void*)0
            pthread_join(thread1_, nullptr);
        }
        thread1_ = 0;
        if (serialDataManager.sp_ != nullptr) {
            serialDataManager.sp_->close();
            delete serialDataManager.sp_;
        }
        serialDataManager.sp_ = nullptr;
    }
    catch(...) {

    }
    printf("quit Ntrip thread\n");
}

int GNSS::ParseRawData(char *message, int len) {
    PosSolver *solver = new PosSolver(svsManager, &rtkManager, this);
    memcpy(solver->raw, message, len);

    solver->CalcuPosition();
//    pthread_create(&threadPos, nullptr, PositionThread, solver);

    return 1;
}

void* GNSS::ThreadAdapterGNSS(void *__sData) {
    auto _sData= ( SerialData* ) __sData;
    _sData->StartCapture(_sData->serialPort_,_sData->baudRate);
    return nullptr;
}

void* GNSS::ThreadAdapterQianXun(void *__rtk) {
    auto _rtk= ( NtripRTK* ) __rtk;
    _rtk->RecvThread();
    return nullptr;
}

void* GNSS::PositionThread(void *__pos) {
    auto _pos = ( PosSolver* ) __pos;
    _pos->CalcuPosition();
}