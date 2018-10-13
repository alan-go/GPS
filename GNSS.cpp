#include "GNSS.h"
#include "EphemSp3.h"

GNSS::GNSS() :useGPS(1),useBeiDou(1),useQianXun(1),isPositioned(false),ephemType(0),logOpen(0){
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

GNSS::~GNSS() {
    for(SerialData* seri:serialDataManager)delete seri;
}

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

    kalmanSolver.InitKalman(this);

}

int GNSS::StartGNSS(const std::string &serial_port, const unsigned int baudRate) {
    //    rtk_protocol_ = rtk_protocol;

}

int GNSS::StartGNSS() {
    sprintf(timeName,"%02d%02d_%02d_%02d",utcTime->tm_mon+1,utcTime->tm_mday,utcTime->tm_hour,utcTime->tm_min);

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
    for(SerialData* seri:serialDataManager){
        pthread_create(&(seri->thread), nullptr, ThreadAdapterSerial, seri);
        sleep(2);
    }

    return 1;
}

int GNSS::StopGNSS() {
    try {
        if (rtkManager.sock_ > 0) {
            close(rtkManager.sock_);
        }
        rtkManager.sock_ = 0;
        rtkManager.stopRTK = true;
        if (thread2_ > 0) {//(void*)0
            pthread_join(thread2_, nullptr);
        }
        thread2_ = 0;
    }
    catch(...) { }
    for (SerialData* seri:serialDataManager)
        seri->StopDevice();

    printf("quit Ntrip thread\n");
}

int GNSS::ParseRawData(char *message, int len) {
    if(len>2048)return -1;
    int sigId =1;
    char raw[2048];
    vector<SV*> svsVisable;
    memcpy(raw,message+6,len-6);

    double rcvtow = *(double*)raw;
    int week = *(uint16_t*)(raw+8);
    int numMeas = *(u_int8_t*)(raw+11);
    GnssTime rTime(week,rcvtow);
    printf("prepare rawdata ,len = %d numMesa=%d\n",len,numMeas);
    if(0==numMeas)return -1;

    for(u_int8_t n = 0;n<numMeas;n++){
        int n32 = n*32;
        uint8_t gnssId = *(uint8_t *)(raw+36+n32);
        uint8_t svId = *(uint8_t *)(raw+37+n32);
        double prMes = *(double*)(raw+16+n32);
        double cpMes = *(double*)(raw+24+n32);
        double doMes = *(float *)(raw+32+n32);
        cpMes *= GetFreq(SysType(gnssId),sigId,1);
        SV* sv = svsManager.GetSv(SysType(gnssId),svId);
        if(sv== nullptr){
            printf("Wrong sv Id! %d,%d\n", gnssId, svId);
            continue;
        }
        sv->AddMmeasure(new Measure(rTime, prMes, cpMes, doMes));
        svsVisable.push_back(sv);
    }

    Test(svsVisable);
    return 0;
}

int GNSS::Test(vector<SV *> svs) {
    printf("coutnt %d\n", count);
    if(++count<100) return -1;
    PosSolver solverSingle(svsManager, &rtkManager, this);

    if(0==solverSingle.PositionSingle(svs)){
        solverSingle.soltion.Show("###Single###");
        AddPosRecord(solverSingle.soltion);
    }

    PosSolver solverRtk(svsManager, &rtkManager, this);
    if(0== solverRtk.PositionRtk(svs))
        solverRtk.soltion.Show("###RTK###");

//    if(0== kalmanSolver.PositionKalman(svs)){
//        solverRtk.soltion.Show("Kalman solution:::");
//    }

}
int GNSS::Peform(vector<SV *> svs) {
    printf("coutnt %d\n", count);
    if(++count<100) return -1;
    PosSolver *solver = new PosSolver(svsManager, &rtkManager, this);
//    if(useQianXun){
////    if(useQianXun&&isPositioned){
//        if(isPositioned)solver->PositionRtkKalman();
////        else solver->PositionRtk();
//        else solver->PositionSingle();
//    } else{
//        solver->PositionSingle();
//    }
    //todo:多线程算电离层会出错
//    pthread_create(&threadPos, nullptr, PositionThread, solver);

}
void* GNSS::ThreadAdapterSerial(void *__sData) {
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
//    _pos->PositionSingle();
}

int GNSS::AddPosRecord(Solution record) {
    records.push_front(record);
    if(records.size()>10240)records.pop_back();
}

