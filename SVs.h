#ifndef SVS_H
#define SVS_H

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;
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
        uint32_t toe,toeF2,toeF3;
        double sqrtA, e,i0,Omega0,M0;
        double Cus,Cuc,Cis,Cic,Crs,Crc,dtn,omega;
        double OmegaDot;
        double IDOT;
        Orbit():toe(0),sqrtA(0),e(0),i0(0),Omega0(0),M0(0), Cus(0),Cuc(0),Cis(0),Cic(0),Crs(0),Crc(0),dtn(0),omega(0),
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
    uint32_t SatH1,URAI;
    uint32_t SOW,WN;
    ionosphere ino;
    uint32_t AODC;//BeiDou
    uint32_t IODC;//GPS
    double toc,a0,a1,a2,a3;
    double TGD,TGD1,TGD2;
    Orbit orbit;
    double prMes, cpMes, doMes;
    double I,T;

    Vector3d position;
    double tsv,tsDelta,tsReal;
public:
    SV();
    ~SV();
    bool JudgeUsable(bool useBeiDou, bool useGps);
    bool CalcuECEF(double rcvtow);
    bool CalcuTime(double rcvtow);
    void PrintInfo(int printType);
    virtual int DecodeSubFrame(uint32_t* dwrds) = 0;
    //head 指32bit()中的头bit（范围：1-32）
    inline uint32_t Read1Word(uint32_t word, int length, int head, bool isInt = false);
    inline uint32_t Read2Word(uint32_t word0,int length0, int head0, uint32_t word1, int length1, int head1, bool isInt = false);
};

class BeiDouSV:public SV{

public:
    int DecodeSubFrame(uint32_t* dwrds){
        if(isBeiDouGEO)DecodeD2(dwrds);
        else DecodeD1(dwrds);
    }
    int DecodeD1(uint32_t* dwrds);
    int DecodeD2(uint32_t* dwrds);

};

class GpsSV:public SV{
public:
    int DecodeSubFrame(uint32_t* dwrds);

};

class SVs{
public:
    GpsSV svGpss[30];
    BeiDouSV svBeiDous[37];
    GNSS* gnss;


public:
    SVs();
    SVs(GNSS* gnss);
    ~SVs();
    void UpdateEphemeris(char * subFrame);

private:
};

#endif