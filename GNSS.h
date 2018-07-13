#ifndef GNSS_H
#define GNSS_H

#include <iostream>
#include <string>
#include <thread>
#include <mutex>
#include <queue>
#include <chrono>
#include <unistd.h>
#include <boost/asio.hpp>

#include <arpa/inet.h>
#include <sys/socket.h>
#include <netinet/in.h>

#include "SVs.h"
#include "PosSolver.h"
#include "NtripRTK.h"
#include "SerialData.h"

class GNSS{
public:
    NtripRTK rtkServer;
    SerialData serialDataManager;
    std::string serialPort_, serverIP_, rtk_protocol_;
    unsigned short port_;
    int sock_;


public:
    GNSS();

    ~GNSS();

    int StartGNSS(const std::string &serial_port, const unsigned int baudRate);

    int StopGNSS();

private:
    pthread_t thread1_, thread2_;
    boost::asio::serial_port *sp_;
    bool stop_rtk_;
private:
    static void *ThreadAdapterGNSS(void *__this);

    static void *ThreadAdapterQianXun(void *__this);

};
#endif