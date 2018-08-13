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
class MSM4Cell{
public:
    double df400,df401;
    double df420,df402;
    double df403;

    MSM4Cell():df400(0),df401(0),df420(0),df402(0),df403(0){}
};
class MSM4data{
public:
    double rtktime;
    double df397, df398;
    vector<int>sigs;
    MSM4Cell sigData[32];

    MSM4data():rtktime(0),df397(0),df398(0){}
};

class NtripRTK{
public:
    GNSS *gnss;
    bool stopRTK;
    std::string serverIP_;
    unsigned short port_;
    int sock_;

    vector<MSM4data*> rtkDataGps[NGPS],rtkDataBeiDou[NBeiDou];
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
    int ParaseMSM4(char *bufferRTK, SV::SvType type);
    bool NtripLogin(const std::string &rtk_protocol);
    void RecvThread();
    int TestParase(char *bufferRecv,int recvLength);
    MSM4data* GetRtkRecord(int satInd,int timeInd, SV::SvType type);

private:
    boost::crc_basic<24> crc24Q;

private:
    int AddRtkRecord(MSM4data* data,SV::SvType type, int id);
    inline uint32_t NetToHost32(char *begin,int head,int length);
    int ParaseRtk32_1005(char * buffer);
};

#endif