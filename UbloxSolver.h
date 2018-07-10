#ifndef UBLOX_SOLVER_H
#define UBLOX_SOLVER_H

#include <iostream>
#include <cmath>
#include <vector>
#include <thread>
#include <Eigen/Core>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;

constexpr static double M_miu = 3.986004418e14;
constexpr static double Omega_e = 7.2921150e-005;
constexpr static double Light_speed = 299792358;
constexpr static double Earth_a = 6378137.0;  //地球长半轴
constexpr static double Earth_f = 3.352810664747481e-003;   //基准椭球体的极扁率  f = 1/298.257223563
constexpr static double Earth_ee = 6.694379990141317e-003;   //偏心率e   e^2 = f(2-f)




class SvInfo
{
public:
    enum SvType{
        GPS = 0,
        BeiDou = 3
    };
    struct Orbit{
        uint32_t toe,toeF2,toeF3;
        double sq_a, e,i0,Omega0,M0;
        double Cus,Cuc,Cis,Cic,Crs,Crc,dtn,omega;
        double OmegaDot;
        double IDOT;
        Orbit():toe(0),sq_a(0),e(0),i0(0),Omega0(0),M0(0),
        Cus(0),Cuc(0),Cis(0),Cic(0),Crs(0),Crc(0),dtn(0),omega(0),
                IDOT(0),OmegaDot(0){}
    };
    struct ionosphere{
        double a0,a1,a2,a3;
        double b0,b1,b2,b3;
        ionosphere():a0(0),a1(0),a2(0),a3(0),b0(0),b1(0),b2(0),b3(0){}
    };

    bool page1OK,page2OK,page3OK,pageOK;
    bool isBeiDouGEO = false;
    bool open = true;
    SvType type;
    int svId;
    uint32_t SOW,WN;
    ionosphere ino;
    uint32_t AODC;
    double toc,a0,a1,a2,a3;
    double TGD1,TGD2;
    uint32_t AODE;
    Orbit orbit;
    double prMes, cpMes, doMes;
    double I,T;

    Vector3d position;
    double ts,tsDelta,tsReal;
public:
    SvInfo();
    ~SvInfo();
    bool CalcuECEF(double rcvtow);
    bool CalcuTime(double rcvtow);
    void PrintInfo(int printType);
};

class UbloxSolver
{
public:
    double rcvtow;
    u_int8_t numMeas = 0;

    Vector4d rxyzt,rxyzOld;//ECEF position of receiver
//    double longitude,latitude,height;
    Vector3d LLA;
    SvInfo GPSSVs[32],GPSSVsCopy[32];
    SvInfo BeiDouSVs[37],BeiDouSVsCopy[37];
    vector<SvInfo*> visibleSvs;
    vector<SvInfo> SvsForCalcu;
    int numtemp = 0;
    bool isCalculating = false;
    bool useGPS = true;
    bool useBeiDou = true;

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
    bool CalcuSvTime();
    //head 指32bit中的头bit（范围：1-32）
    inline uint32_t Read1Word(uint32_t word, int length, int head = 2,bool isInt = false);
    inline uint32_t Read2Word(uint32_t* word,int length0, int head0,int length1, int head1 = 2,bool isInt = false);
    inline uint32_t Read3Word(uint32_t* word, int length0, int head0,int length1,int head1, int length2, int head2, bool isInt = false);


public:
    UbloxSolver();
    ~UbloxSolver();
    bool ParseRawData(char* message);
    bool ParseBstSubFrame(char* message);
    bool XYZ2LLA(Vector3d XYZ,Vector3d &LLA);
};

#endif //UBLOX_SOLVER_H
