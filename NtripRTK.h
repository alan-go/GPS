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
#include <mutex>


#include "SVs.h"

class GNSS;
class MSM4Cell{//每个星的每一个信号有一个MSM4Cell
public:
    double df400,df401;
    double df420,df402;
    double df403;
    double prMes,cpMes;
    MSM4Cell():df400(0),df401(0),df420(0),df402(0),df403(0),prMes(0),cpMes(0){}
};

class MSM4data{//每个星有一个MSM4data
public:
    double rtktime;
    uint32_t refID;
    double df397, df398, prRough;
    vector<int>sigs;
    MSM4Cell sigData[32];
    bool isOk;
    MSM4data():rtktime(0),refID(0),df397(0),df398(0),prRough(0),isOk(0){}
};

class NtripRTK{
public:
    GNSS *gnss;
    bool stopRTK;
    std::string serverIP_;
    unsigned short port_;
    int sock_;
    //vector<MSM4data*>里面是不同时刻的数据，0是最近的记录
//    vector<MSM4data*> rtkDataGps[NGPS],rtkDataBeiDou[NBeiDou];
    uint32_t refStationId;
    bool isPhysicalStation;
    bool singleReceiver;
    uint32_t ITRFyear;
    bool supportGPS,supportBeiDou,supportGLONASS,supportGalileo;
    uint32_t quaterCycle;
    Eigen::Vector3d ECEF_XYZ;
    char ggaDefault[128];
    std::mutex mtxData;

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