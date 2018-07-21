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




#include "SVs.h"
#include "PosSolver.h"
#include "NtripRTK.h"
#include "SerialData.h"

class GNSS{
public:
    Eigen::Matrix<double,3,1> xyz, xyzOld;
    double tu, tuBeiDou, tuGps;
    SVs svsManager;
    SerialData serialDataManager;
    NtripRTK rtkManager;
    std::string serialPort_, rtk_protocol_;

    bool useBeiDou,useGPS;

public:
    GNSS();

    ~GNSS();

    int StartGNSS(const std::string &serial_port, const unsigned int baudRate);

    int StartGNSS(const std::string &fileName);

    int StopGNSS();

    int ParseRawData(char * message, int len);

private:
    pthread_t thread1_, thread2_, threadPos;
private:
    static void *ThreadAdapterGNSS(void *__sData);

    static void *ThreadAdapterQianXun(void *__rtk);

    static void *PositionThread(void *__pos);

};
#endif