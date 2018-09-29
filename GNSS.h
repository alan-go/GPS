#ifndef GNSS_H
#define GNSS_H

#include "CommonInclude.h"
#include "SVs.h"
#include "PosSolver.h"
#include "NtripRTK.h"
#include "SerialData.h"
#include "GnssTime.h"

class GNSS{
public:
    struct PosRcd{
        double rcvtow;
        Eigen::Vector3d xyz, vxyz, lla;
        PosRcd(double tow,Vector3d xyz,Vector3d vxyz, Vector3d lla):rcvtow(tow),xyz(xyz),vxyz(vxyz),lla(lla){}
    };
    Eigen::Vector3d xyzDefault, llaDefault;
//    Eigen::Vector3d LLA, LLAOld;
//    vector<PosRcd> records;
    deque<PosRcd,Eigen::aligned_allocator<Eigen::Vector3d>> records;
    double cycle[Nsys-1][Nxxs],PB[Nsys-1][Nxxs],sigmaCy[Nsys-1][Nxxs],sigmaPr[Nsys-1][Nxxs];
    Matrix<double,6,6> Pxv;
//    Matrix<double,Ngps,1> cycleGPS,sigmaGPSCy,sigmaGPSPr;
    double tu, tuBeiDou, tuGps;
    SVs *svsManager;
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
    std::FILE *log,*logRTK;
    struct tm *utcTime;

public:
    GNSS();

    ~GNSS();

    int Init(int ephem,bool qianXun, bool bds, bool gps);

    int StartGNSS(const std::string &serial_port, const unsigned int baudRate);

    int StartGNSS(const std::string &fileName);

    int StopGNSS();

    int ParseRawData(char * message, int len);

    int AddPosRecord(PosRcd record);

private:
    pthread_t thread1_, thread2_, threadPos;
private:
    static void *ThreadAdapterGNSS(void *__sData);

    static void *ThreadAdapterQianXun(void *__rtk);

    static void *PositionThread(void *__pos);

};
#endif