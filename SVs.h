#ifndef SVS_H
#define SVS_H

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;
#define NGPS 64
#define NBeiDou 64
//PI ???GPs_PI
constexpr static double GPS_PI = 3.1415926535898;
constexpr static double M_miu = 3.986004418e14;
constexpr static double Omega_e = 7.2921150e-005;//GPS?Beidou
constexpr static double Light_speed = 299792358;
constexpr static double Earth_a = 6378137.0;  //地球长半轴
constexpr static double Earth_f = 3.352810664747481e-003;   //基准椭球体的极扁率  f = 1/298.257223563
constexpr static double Earth_ee = 6.694379990141317e-003;   //偏心率e   e^2 = f(2-f)

class GNSS;
class SV{
public:
    enum SvType{
        GPS = 0,
        BeiDou = 3
    };
    struct Orbit{
        uint32_t AODE;//BeiDou
        uint32_t IODE;//GPS
        uint32_t toe;
        double sqrtA, e,i0,Omega0,M0;
        double Cus,Cuc,Cis,Cic,Crs,Crc,dtn,omega;
        double OmegaDot;
        double IDOT;

        uint32_t CucHigh,CucLow;
        uint32_t CicHigh,CicLow;
        uint32_t eHigh,eLow;
        uint32_t toeHigh,toeLow;
        uint32_t i0High,i0Low;
        uint32_t OmegaDotHigh,OmegaDotLow;
        uint32_t omegaHigh,omegaLow;
        Orbit():toe(0),sqrtA(0),e(0),i0(0),Omega0(0),M0(0), Cus(0),Cuc(0),Cis(0),Cic(0),Crs(0),Crc(0),dtn(0),omega(0),
                IDOT(0),OmegaDot(0){}
    };
    struct ionosphere{
        double a0,a1,a2,a3;
        double b0,b1,b2,b3;
        ionosphere():a0(0),a1(0),a2(0),a3(0),b0(0),b1(0),b2(0),b3(0){}
    };
    struct SignalData{
        double df400,df401;
        uint32_t df402,df420,df403;
        SignalData():df400(0),df401(0),df402(0),df420(0),df403(0){}
    };

    int8_t bstEphemOK[10];

    bool isBeiDouGEO;
    bool open, measureGood, elevGood;
        SvType type;
    int svId;
    uint32_t SatH1,URAI;
    uint32_t SOW,WN;
    ionosphere ino;
    uint32_t AODC;//BeiDou
    uint32_t IODC;//GPS
    double toc,a0,a1,a2,a3;
    uint32_t a1High,a1Low;
    double TGD1,TGD2;
    Orbit orbit;
    double prMes, cpMes, doMes;
    double I,T;

    Vector3d position,sLLA;
    double tsv,tsDelta,tsReal;

    double elevationAngle,azimuthAngle;
    //for RTK:
    double df397,df398;
public:
    SV();
    ~SV();
    bool JudgeUsable(bool useBeiDou, bool useGps);
    bool MeasureGood();
    bool ElevGood();
    bool CalcuECEF(double rcvtow);
    bool CalcuTime(double rcvtow);
    int CalcuelEvationAzimuth(Vector3d receiverPosition, Vector3d LLA);
    int CalcuTroposhphere(double elev,double azim);
    int CalcuInoshphere(double elev,double azim,Vector3d LLA,double time);
    int CorrectIT(Vector3d receiverPosition, Vector3d LLA,double time);
    void PrintInfo(int printType);
//    virtual int DecodeSubFrame(uint32_t* dwrds) = 0;
    virtual SignalData* SignalTable(int index) = 0;

    //head 指32bit()中的头bit（范围：1-32）
    inline uint32_t Read1Word(uint32_t word, int length, int head, bool isInt = false);
    inline uint32_t Read2Word(uint32_t word0,int length0, int head0,
            uint32_t word1, int length1, int head1, bool isInt = false);
    inline uint32_t Read3Word(uint32_t word0,int length0, int head0,
            uint32_t word1, int length1, int head1, uint32_t word2, int length2, int head2, bool isInt = false);
};

class BeiDouSV:public SV{
public:
    //signal number:2,3,4;  8,9,10;  14,15,16;
    SignalData B1_2I,B1_2Q,B1_2X,B3_6I,B3_6Q,B3_6X,B2_7I,B2_7Q,B2_7X;
public:
    BeiDouSV();
    ~BeiDouSV();
    int DecodeSubFrame(uint32_t* dwrds){
        if(isBeiDouGEO)DecodeD2Frame1(dwrds);
        else DecodeD1(dwrds);
        return 1;
    }
    int DecodeD1(uint32_t* dwrds);
    int DecodeD2Frame1(uint32_t *dwrds);
    SignalData* SignalTable(int index);

};

class GpsSV:public SV{
public:
    SignalData L1_1C,L1_1P,L1_1W,L2_2C,L2_2P,L2_2W,L2_2S,L2_2L,L2_2X,L5_5I,L5_5Q,L5_5X;
public:
    GpsSV();
    ~GpsSV();
    int DecodeSubFrame(uint32_t* dwrds);
    SignalData* SignalTable(int index);
};

class SVs{
public:
    GpsSV svGpss[NGPS];
    BeiDouSV svBeiDous[NBeiDou];
    GNSS* gnss;


public:
    SVs();
//    SVs(GNSS* gnss);
    ~SVs();
    void UpdateEphemeris(char * subFrame);
    SV* SatTable(SV::SvType type,int ind);
private:
};

#endif