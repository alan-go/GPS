#include "GNSS.h"
#include "EphemSp3.h"

GNSS::GNSS() :tu(0),tuBds(0),tuGps(0),useGPS(1),useBeiDou(1),useQianXun(1),isPositioned(false),ephemType(0),logOpen(0){
    serialDataManager.gnss = this;
    rtkManager.gnss = this;

    rtkManager.serverIP_ = "60.205.8.49";
    rtkManager.port_ = 8002;
    xyzDefault<<-2.17166e+06, 4.38439e+06, 4.07802e+06;
    llaDefault<<40.0*R2D, 116.345*R2D, 59.07;

    time_t timeNow = time(NULL);
    utcTime=gmtime(&timeNow);
    memset(svMaskBds,1,Nbds * sizeof(int));
    memset(svMaskGps,1,Ngps * sizeof(int));
}

GNSS::~GNSS() {}

int GNSS::Init(int ephem, bool qianXun, bool bds, bool gps) {
    useGPS = gps;
    useBeiDou = bds;
    svsManager.gnss = this;
    svsManager.InitAlloc();
    svsManager.SetOpen(bds,gps);
    ephemType = ephem;
    useQianXun = qianXun;
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
    //todo
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
    if(len>256)return -1;
    char raw[256];
    vector<SV*> svsVisable;
    memcpy(raw,message+6,len-6);

    double rcvtow = *(double*)raw;
    int week = *(uint16_t*)(raw+8);
    int numMeas = *(u_int8_t*)(raw+11);
    GnssTime rTime(week,rcvtow);
    printf("prepare rawdata , numMesa=%d\n",numMeas);
    if(0==numMeas)return -1;

    for(u_int8_t n = 0;n<numMeas;n++){
        int n32 = n*32;
        uint8_t gnssId = *(uint8_t *)(raw+36+n32);
        uint8_t svId = *(uint8_t *)(raw+37+n32);
        double prMes = *(double*)(raw+16+n32);
        double cpMes = *(double*)(raw+24+n32);
        double doMes = *(float *)(raw+32+n32);
        SV* sv = svsManager.GetSv(SysType(gnssId),svId);
        if(sv== nullptr){
            printf("Wrong sv Id! %d,%d\n", gnssId, svId);
            continue;
        }
        if((rTime-sv->measureDat.front()->time)>1.5)sv->trackCount=0;
        else sv->trackCount++;
        Measure *mesur = new Measure(rTime,prMes,cpMes,doMes);
        sv->measureDat.push_front(mesur);
        svsVisable.push_back(sv);
    }
    Test(svsVisable);
    return 0;
}

int GNSS::Test(vector<SV *> svs) {
    printf("coutnt %d\n", count);
    if(++count<100) return -1;
    PosSolver solverSingle(svsManager, &rtkManager, this);
    PosSolver solverSingle2sys(svsManager, &rtkManager, this);
    PosSolver solverRtk(svsManager, &rtkManager, this);
    PosSolver solverKalman(svsManager, &rtkManager, this);

//    solverSingle.
    if(useQianXun){
//    if(useQianXun&&isPositioned){
        if(isPositioned)solverSingle.PositionRtkKalman();
//        else solverSingle->PositionRtk();
        else solverSingle.PositionSingle();
    } else{
        solverSingle.PositionSingle();
    }

}
int GNSS::Peform(vector<SV *> svs) {
    printf("coutnt %d\n", count);
    if(++count<100) return -1;
    PosSolver *solver = new PosSolver(svsManager, &rtkManager, this);
    if(useQianXun){
//    if(useQianXun&&isPositioned){
        if(isPositioned)solver->PositionRtkKalman();
//        else solver->PositionRtk();
        else solver->PositionSingle();
    } else{
        solver->PositionSingle();
    }
    //todo:多线程算电离层会出错
//    pthread_create(&threadPos, nullptr, PositionThread, solver);

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

int GNSS::AddPosRecord(Solution record) {
    records.push_front(record);
    if(records.size()>10240)records.pop_back();
}

