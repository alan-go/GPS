#ifndef GNSS_H
#define GNSS_H

#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <mutex>
#include <queue>
#include <chrono>
#include <unistd.h>
#include <stdio.h>



#include "SVs.h"
#include "PosSolver.h"
#include "NtripRTK.h"
#include "SerialData.h"

class GNSS{
    struct PosRcd{
        double rcvtow;
        Eigen::Vector3d xyz,lla;
        PosRcd(double tow,Vector3d xyz, Vector3d lla):rcvtow(tow),xyz(xyz),lla(lla){}
    };
public:
//    Eigen::Vector3d xyz, xyzOld;
//    Eigen::Vector3d LLA, LLAOld;
    vector<PosRcd> records;
    double tu, tuBeiDou, tuGps;
    SVs svsManager;
    SerialData serialDataManager;
    NtripRTK rtkManager;
    std::string serialPort_, rtk_protocol_;

    bool useBeiDou,useGPS,useQianXun;
    bool isPositioned;


public:
    GNSS();

    ~GNSS();

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