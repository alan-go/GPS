#ifndef GNSS_H
#define GNSS_H

#include "CommonInclude.h"
#include "Svs.h"
#include "PosSolver.h"
#include "NtripRTK.h"
#include "SerialData.h"
#include "GnssTime.h"


class GNSS{
public:
    bool stop{false};
    Eigen::Vector3d xyzDefault, llaDefault;
    Eigen::Vector3d xyz00,xyzUBX,xyzRAC;
    deque<Solution,Eigen::aligned_allocator<Eigen::Vector3d>> records;
    deque<Solution,Eigen::aligned_allocator<Eigen::Vector3d>> solSingleNew,solSingles,solRTKs,solKalmans,solKalDops,solSigKals;
    Solution solSingle,solRTK,solKalman,solUBX,solRAC;
    double cycle[Nsys-1][Nxxs],PB[Nsys-1][Nxxs],sigmaCy[Nsys-1][Nxxs],sigmaPr[Nsys-1][Nxxs];
    Matrix<double,6,6> Pxv;
//    Matrix<double,Ngps,1> cycleGPS,sigmaGPSCy,sigmaGPSPr;
    SvAll svsManager;
    vector<SerialData*> serialDataManager;
    NtripRTK rtkManager;
    std::string rtk_protocol_;

    bool useBeiDou,useGPS,useQianXun;
    int svMaskBds[Nbds], svMaskGps[Ngps];
    int ephemType;//0:broadcast,1:sp3
//    bool isPositioned;
    int count{0};
//    char logFile[64];
    bool logOpen{1};
    std::FILE *log,*logRTK,*logDebug,*logTu;
    std::FILE *logSingle,*logSingleNew,*logKalman;
    struct tm *utcTime;
    char timeName[128];
    PosSolver kalmanSolver,solver,kalDoppSolv;
    GnssTime rTime;

public:
    GNSS();

    ~GNSS();

    int Init(int ephem,bool qianXun, bool bds, bool gps);

    int StartGNSS(const std::string &serial_port, const unsigned int baudRate);

    int StartGNSS();

    int StopGNSS();

    int ParseRawData(char * message, int len);

    int AddPosRecord(Solution record);

    int AddSerial(int id,int type,string name, unsigned int baudRate,bool logOpen,bool paraseData)
        {serialDataManager.push_back(new SerialData(id, type,name,baudRate,logOpen,paraseData,this));};

    SerialData* GetSerial(int id){
        for(SerialData* seri:serialDataManager)
            if(id==seri->id)
                return seri;
        return nullptr;
    }


    int Peform(vector<SV*> svs);
    int Test(vector<SV*> svs);

private:
    pthread_t thread1_, thread2_, threadPos, threadSp3_;
private:
    static void *ThreadAdapterSerial(void *__sData);

    static void *ThreadAdapterQianXun(void *__rtk);

    static void *PositionThread(void *__pos);

    static void *ThreadAdapterSp3(void *__gnss);

};
#endif