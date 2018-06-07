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
        uint32_t toe,toeF2,toeF3;
        double sq_a, e,i0,Omega0,w,M0;
        double Cus,Cuc,Cis,Cic,Crs,Crc,dtn,omega,I;
        double OmegaDot;
        double IDOT;
        Orbit():toe(0),sq_a(0),e(0),i0(0),Omega0(0),w(0),M0(0),
        Cus(0),Cuc(0),Cis(0),Cic(0),Crs(0),Crc(0),dtn(0),omega(0),
        I(0),IDOT(0),OmegaDot(0){}
    };
    struct ionosphere{
        double a0,a1,a2,a3;
        double b0,b1,b2,b3;
        ionosphere():a0(0),a1(0),a2(0),a3(0),b0(0),b1(0),b2(0),b3(0){}
    };
public:
    SvId id;
    uint32_t SOW,WN;
    ionosphere ino;
    uint32_t AODC;
    double toc,a0,a1,a2,a3;
    double TGD1,TGD2;
    uint32_t AODE;
    Orbit orbit;
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
    inline uint32_t Read1Word(uint32_t word, int length, int head = 2,bool isInt = false);
    inline uint32_t Read2Word(uint32_t* word,int length0, int head0,int length1, int head1 = 2,bool isInt = false);
    inline uint32_t Read3Word(uint32_t* word, int length0, int head0,int length1,int head1, int length2, int head2, bool isInt = false);


public:
    UbloxSolver();
    ~UbloxSolver();
    bool ParseRawData(char* message);
    bool ParseBstSubFrame(char* message);
};

#endif //UBLOX_SOLVER_H
