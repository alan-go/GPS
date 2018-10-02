#include "GNSS.h"
#include "EphemSp3.h"

GNSS::GNSS() :tu(0),tuBeiDou(0),tuGps(0),useGPS(1),useBeiDou(1),useQianXun(1),isPositioned(false),ephemType(0),logOpen(0){
    serialDataManager.gnss = this;
    rtkManager.gnss = this;

    rtkManager.serverIP_ = "60.205.8.49";
    rtkManager.port_ = 8002;
    xyzDefault<<-2.17166e+06, 4.38439e+06, 4.07802e+06;
    llaDefault<<40.0*GPS_PI/180.0, 116.345*GPS_PI/180.0, 59.07;

    time_t timeNow = time(NULL);
    utcTime=gmtime(&timeNow);
    memset(svMaskBds,1,Nbds * sizeof(int));
    memset(svMaskGps,1,Ngps * sizeof(int));
}

GNSS::~GNSS() {}

int GNSS::Init(int ephem, bool qianXun, bool bds, bool gps) {
    svsManager = new SVs(bds,gps);
    svsManager->gnss = this;
    ephemType = ephem;
    useQianXun = qianXun;
    useGPS = gps;
    useBeiDou = bds;
    if(1==ephemType){
        if(0==EphemSp3::ReadSp3File("/home/alan/Desktop/hour20175_19.sp3",svsManager))return 0;
    } else{
        printf("ReadSp3File Failed. \n");
//        return -1;
    }

    for (int i = 0; i < Nsys - 1; ++i) {
        for (int j = 0; j < Nxxs; ++j) {
//            cycle[i][j]=40.0;
            cycle[i][j]=10;
            sigmaCy[i][j]=0.15;
            sigmaPr[i][j]=1;
            PB[i][j]=100;
        }
    }
    Pxv = 10000*MatrixXd::Identity(6,6);
}

int GNSS::StartGNSS(const std::string &serial_port, const unsigned int baudRate) {
    //    rtk_protocol_ = rtk_protocol;
    serialDataManager.serialPort_ = serial_port;
    serialDataManager.baudRate = baudRate;
    serialDataManager.stopCapture = false;
    sprintf(serialDataManager.saveName,"../data/%02d%02d_%02d_%02d.data",
            utcTime->tm_mon+1,utcTime->tm_mday,utcTime->tm_hour,utcTime->tm_min);
    sprintf(rtkManager.saveName,"../data/%02d%02d_%02d_%02d.rtk",
            utcTime->tm_mon+1,utcTime->tm_mday,utcTime->tm_hour,utcTime->tm_min);

    if(useQianXun){
        if(!rtkManager.NtripLogin(rtk_protocol_)) {
            printf("cannot login ntrip server\n");
            exit(0);
        }
        rtkManager.SentGGA(rtkManager.ggaDefault,strlen(rtkManager.ggaDefault));

        pthread_create(&thread2_, nullptr, ThreadAdapterQianXun, &rtkManager);
        sleep(2);
    }

    //todo : for temmp
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
    printf("coutnt %d\n", count);
    if(++count<1500)
        return -1;
    PosSolver *solver = new PosSolver(svsManager, &rtkManager, this);
    memcpy(solver->raw, message, len);
    if(useQianXun){
//    if(useQianXun&&isPositioned){
        if(isPositioned)solver->PositionRtkKalman();
        else solver->PositionRtk();
    } else{
        solver->PositionSingle();
    }
    //todo:多线程算电离层会出错
//    pthread_create(&threadPos, nullptr, PositionThread, solver);

    return 1;
}

void* GNSS::ThreadAdapterGNSS(void *__sData) {
    auto _sData= ( SerialData* ) __sData;
    _sData->StartCapture(_sData->serialPort_,_sData->baudRate,_sData->saveName);
    return nullptr;
}

void* GNSS::ThreadAdapterQianXun(void *__rtk) {
    auto _rtk= ( NtripRTK* ) __rtk;
    _rtk->RecvThread();
    return nullptr;
}

void* GNSS::PositionThread(void *__pos) {
    auto _pos = ( PosSolver* ) __pos;
    _pos->PositionSingle();
}

int GNSS::AddPosRecord(GNSS::PosRcd record) {
    records.push_front(record);
    if(records.size()>10240)records.pop_back();
}

