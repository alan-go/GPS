#ifndef NTRIPRTK_H
#define NTRIPRTK_H

#include <iostream>
#include <string>
#include <unistd.h>

#include <arpa/inet.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <boost/asio.hpp>
#include <boost/crc.hpp>


#include "SVs.h"

class GNSS;
class NtripRTK{
public:
    GNSS *gnss;
    bool stopRTK;
    std::string serverIP_;
    unsigned short port_;
    int sock_;

    int satGps[NGPS],satBeidou[NBeiDou];
    int sigGps[NGPS],sigBeidou[NBeiDou];
    uint32_t refStationId;
    bool isPhysicalStation;
    bool singleReceiver;
    uint32_t ITRFyear;
    bool supportGPS,supportBeiDou,supportGLONASS,supportGalileo;
    uint32_t quaterCycle;
    Eigen::Vector3d ECEF_XYZ;

public:
    NtripRTK();

    ~NtripRTK();

    void UpdateGGA();
    int SentGGA(const char *bufferGGA, int length);
    int ParaseMSM4(char *bufferRTK, SV *sv, int *sat, int *sig);
    bool NtripLogin(const std::string &rtk_protocol);
    void RecvThread();
    int TestParase(char *bufferRecv,int recvLength);
private:
    boost::crc_basic<24> crc24Q;

private:
    inline uint32_t NetToHost32(char *begin,int head,int length);
    int ParaseRtk32_1005(char * buffer);
};

#endif