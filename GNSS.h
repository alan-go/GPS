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

#include <arpa/inet.h>
#include <sys/socket.h>
#include <netinet/in.h>


#include "SVs.h"
#include "PosSolver.h"
#include "NtripRTK.h"
#include "SerialData.h"

class GNSS{
public:
    Eigen::Matrix<double,3,1> xyz, xyzOld;
    double tu, tuBeiDou, tuGps;
    NtripRTK rtkServer;
    SerialData serialDataManager;
    SVs svsManager;
    NtripRTK rtkManager;
    std::string serialPort_, serverIP_, rtk_protocol_;
    unsigned short port_;
    int sock_;
    bool useBeiDou,useGPS;

public:
    GNSS();

    ~GNSS();

    int StartGNSS(const std::string &serial_port, const unsigned int baudRate);

    int StartGNSS(const std::string &fileName);

    int StopGNSS();

    int ParseRawData(char * message);

private:
    pthread_t thread1_, thread2_, threadPos;
    boost::asio::serial_port *sp_;
    bool stop_rtk_;
private:
    static void *ThreadAdapterGNSS(void *__this);

    static void *ThreadAdapterQianXun(void *__this);

    static void *PositionThread(void *__pos);

};
#endif