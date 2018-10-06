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

    Eigen::Vector3d xyzDefault, llaDefault;
    Eigen::Vector3d xyz00;
//    Eigen::Vector3d LLA, LLAOld;
//    vector<Solution> records;
    deque<Solution,Eigen::aligned_allocator<Eigen::Vector3d>> records;
    double cycle[Nsys-1][Nxxs],PB[Nsys-1][Nxxs],sigmaCy[Nsys-1][Nxxs],sigmaPr[Nsys-1][Nxxs];
    Matrix<double,6,6> Pxv;
//    Matrix<double,Ngps,1> cycleGPS,sigmaGPSCy,sigmaGPSPr;
    double tu, tuBds, tuGps;
    SvAll svsManager;
    SerialData serialDataManager;
    NtripRTK rtkManager;
    std::string serialPort_, rtk_protocol_;

    bool useBeiDou,useGPS,useQianXun;
    int svMaskBds[Nbds], svMaskGps[Ngps];
    int ephemType;//0:broadcast,1:sp3
    bool isPositioned;
    int count{0};
//    char logFile[64];
    bool logOpen;
    std::FILE *log,*logRTK,*logDebug;
    struct tm *utcTime;

public:
    GNSS();

    ~GNSS();

    int Init(int ephem,bool qianXun, bool bds, bool gps);

    int StartGNSS(const std::string &serial_port, const unsigned int baudRate);

    int StartGNSS(const std::string &fileName);

    int StopGNSS();

    int ParseRawData(char * message, int len);

    int AddPosRecord(Solution record);

    int Peform(vector<SV*> svs);
    int Test(vector<SV*> svs);

private:
    pthread_t thread1_, thread2_, threadPos;
private:
    static void *ThreadAdapterGNSS(void *__sData);

    static void *ThreadAdapterQianXun(void *__rtk);

    static void *PositionThread(void *__pos);

};
#endif