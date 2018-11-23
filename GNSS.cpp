#include "GNSS.h"
#include "EphemSp3.h"

GNSS::GNSS() :useGPS(1),useBeiDou(1),useQianXun(1),ephemType(0),logOpen(0){
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
    solver.Init(this);

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
    int tracStat = *(uint8_t *)(raw+12);
    fprintf(logTu,"%f,%d\n",rcvtow,tracStat);
    rTime = GnssTime(week,rcvtow);
    printf("prepare rawdata ,len = %d numMesa=%d\n",len,numMeas);
    if(0==numMeas)return -1;

    for(u_int8_t n = 0;n<numMeas;n++){
        Measure *mes = new Measure();
        mes->time=rTime;
        int n32 = n*32;
        uint8_t gnssId = *(uint8_t *)(raw+36+n32);
        uint8_t svId = *(uint8_t *)(raw+37+n32);
        SV* sv = svsManager.GetSv(SysType(gnssId),svId);
        if(sv== nullptr){
            printf("Wrong sv Id! %d,%d\n", gnssId, svId);
            continue;
        }

        double lambda = GetFreq(SysType(gnssId),sigId,1);
//        double freq = GetFreq(SysType(gnssId),sigId,0);

        mes->prMes = *(double*)(raw+16+n32);
        mes->cpMes = *(double*)(raw+24+n32)*lambda;
        mes->doMes = -(*(float *)(raw+32+n32))*lambda;

        mes->lockTime = *(uint16_t*)(raw+40+n32);
        mes->cno = *(uint8_t *)(raw+42+n32);
        mes->stdevPr = 0.01 *pow(2, *(uint8_t *)(raw+43+n32));
        mes->stdevCp = 0.004*lambda*(*(uint8_t *)(raw+44+n32));
        mes->stdevDo = 0.002*pow(2, *(uint8_t *)(raw+45+n32))*lambda;
        mes->trkStat = *(uint8_t *)(raw+46+n32);

        sv->AddMmeasure(mes);
        svsVisable.push_back(sv);
    }

    Test(svsVisable);
    return 0;
}

int GNSS::Test(vector<SV *> svs) {
    printf("coutnt %d\n", ++count);
    if(count<1000) return -1;
    for(SV*sv:svs){
        sv->FPrintInfo(0);

//        sv->SmoothKalman0();
//        sv->SmoothKalman(1);
//        sv->InterpMeasere(25,2);
    }
//    return 0;

    double tod = fmod(rTime.tow,86400.0);
    double tu_s=solver.solSingle.tu[SYS_GPS]/Light_speed;
//    fprintf(logTu,"%f,%.10f\n",rTime.tow,tu_s*Light_speed);
    solRAC = FindSol(GetSerial(1)->solRaw,tod-tu_s,0.2,"tod");
    xyzRAC=solRAC.xyz;
    solUBX = FindSol(GetSerial(0)->solRaw,tod-tu_s,1,"tod");
    solUBX.Show("###UBX###");
    (solUBX-solRAC).Show("###UBX-RAC###",1);
    solRAC.Show("###RAC###");
    printf("tow,tod insolution= %f,%f\n", rTime.tow,tod-tu_s);
    if(solver.PrepareSVsData(svs))return -1;
//    kalmanSolver.PrepareSVsData(svs);
    printf("mamamaha\n");
    printf("mamamaha\n");
    printf("mamamaha\n");

    if(0==solver.PositionSingle(svs)){
        solver.solSingle.Show("###Single###");
        (solver.solSingle-solRAC).Show("###SIG-RAC###",1);
        AddPosRecord(solver.solSingle);
        solSingles.push_front(solver.solSingle);
//        solver.solSingle.printXYZ2log(logSingle);
    } else{ printf("PositionSingle failed %d\n",count );}

    if(1){
        solver.PosKalSng(svs);
        solver.solKalSigle.Show("###PosKalSng###");
        (solver.solKalSigle-solRAC).Show("###KalSIG-RAC###",1);
        solSigKals.push_front(solver.solSingle);
    }
    return 0;


    if(0==solver.PositionSingleNew(svs)){
        solver.solSingle.Show("###SingleNew###");
        (solver.solSingle-solRAC).Show("###SIGNew-RAC###",1);
//        solver.solSingle.printXYZ2log(logSingleNew);
        solSingleNew.push_front(solver.solSingle);
    }

//    PosSolver solverRtk(svsManager, &rtkManager, this);
    if(0== solver.PositionRtk(svs)){
        solver.solSingle.Show("###RTK###");
        (solver.solSingle-solRAC).Show("###RTK-RAC###",1);
        solRTKs.push_front(solver.solSingle);
    }

//    return 0;
    if(0== kalmanSolver.PositionKalman(svs)){
        kalmanSolver.solSingle.Show("###Kalman###");
        (solver.solSingle-solRAC).Show("###KAL-RAC###",1);
        solKalmans.push_front(kalmanSolver.solSingle);
    }

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

