#ifndef UBLOX_SOLVER_H
#define UBLOX_SOLVER_H

#include <iostream>
#include <cmath>
#include <vector>
#include <thread>
#include <mutex>

using namespace std;

class SvInfo
{
    enum SvId{};
    struct Orbit{
        double sq_a, e,i0,omega0,w,M0;
        double Cus,Cuc,Cis,Cic,Crs,Crc,dtn,Omega,I;
        double toe,idoe;
        Orbit():sq_a(0),e(0),i0(0),omega0(0),w(0),M0(0),
        Cus(0),Cuc(0),Cis(0),Cic(0),Crs(0),Crc(0),dtn(0),Omega(0),
        I(0),toe(0),idoe(0){}
    };
public:
    SvId id;
    Orbit orbit;
    double toc,a0,a1,a2,a3;
    double x,y,z;
    double prMes, cpMes, doMes;
    mutex mtx;
public:
    bool CalcuECEF();
};

class UbloxSolver
{
public:
    double rcvtow;
    u_int8_t numMeas = 0;
    double rx,ry,rz;//ECEF position of receiver
    double longitude,latitude,height;
    SvInfo GPSSVs[32];
    SvInfo BeiDouSVs[37];
    vector<int> GPSIndex,BeiDouIndex;
    int numtemp = 0;
    mutex positionMtx;

private:
    static void* LaunchPositionThread ( void* __this ) {
        auto _this= ( UbloxSolver* ) __this;
        _this->solvePosition();
        return nullptr;
    }

    bool DecodeGpsBroadcast(uint32_t* dwrds,SvInfo* sv);
    bool DecodeBeiDouBroadcastD1(uint32_t* dwrds,SvInfo* sv);
    bool DecodeBeiDouBroadcastD2(uint32_t* dwrds,SvInfo* sv);
    bool solvePosition();
    //head 指32bit中的头bit（范围：1-32）
    inline uint32_t Read1Word(uint32_t word,int head,int length);
    inline uint64_t Read2Word(uint32_t word[2],int head[2],int length[2]);
    inline uint64_t Read3Word(uint32_t word[3],int head[3],int length[3]);

public:
    UbloxSolver();
    ~UbloxSolver();
    bool ParseRawData(char* message);
    bool ParseBstSubFrame(char* message);
};

#endif //UBLOX_SOLVER_H
