#ifndef SVS_H
#define SVS_H

#include "CommonInclude.h"
#include "GnssTime.h"

using namespace std;
using namespace Eigen;


class GNSS;
class MSM4data;
class EphemSp3;
class SV{
public:
    struct Orbit{
        uint32_t AODE;//SYS_BDS
        uint32_t IODE;//SYS_GPS
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

    EphemSp3 *ephemSp3;
    int trackCount;
    bool temp = 1;
    int8_t bstEphemOK[10];
    bool isBeiDouGEO;
    bool open, measureGood, elevGood;
    SysType type;
    int svId;
    uint32_t SatH1,URAI;
    uint32_t SOW,WN;
    ionosphere ino;
    uint32_t AODC;//SYS_BDS
    uint32_t IODC;//SYS_GPS
    double toc,a0,a1,a2,a3;
    uint32_t a1High,a1Low;
    double TGD1,TGD2;
    Orbit orbit;
    //todo:remove
    double prMes, cpMes, doMes;
    double I,T;

    Vector3d position,sLLA;
    double tsv,tsDelta,tsReal;

    double elevationAngle,azimuthAngle;
    //for RTK:
    //vector<MSM4data*>里面是不同时刻的数据，0 always是最近的记录
    deque<MSM4data*> rtkData;
    deque<Measure*> measureDat;
    double prInterp[32], cpInterp[32];
    bool KalmanFirst{1};
public:
    SV();
    SV(int id);
    ~SV();
    bool IsEphemOK(int ephemType, GnssTime time);
    bool MeasureGood();
    bool ElevGood();
    bool IsMaskOn();
    bool CalcuECEF(double tow);
    bool CalcuTime(double tow);
    int CalcuelEvationAzimuth(Vector3d receiverPosition, Vector3d LLA);
    int CalcuTroposhphere(double elev,double azim);
    int CalcuInoshphere(double elev,double azim,Vector3d LLA,double time);
    int CorrectIT(Vector3d pos,double time);
    void PrintInfo(int printType);
    virtual int DecodeSubFrame(uint32_t* dwrds) = 0;
    double InterpRtkData(double time, int sigInd);

    //head 指32bit()中的头bit（范围：1-32）
    inline uint32_t Read1Word(uint32_t word, int length, int head, bool isInt = false);
    inline uint32_t Read2Word(uint32_t word0,int length0, int head0,
            uint32_t word1, int length1, int head1, bool isInt = false);
    inline uint32_t Read3Word(uint32_t word0,int length0, int head0,
            uint32_t word1, int length1, int head1, uint32_t word2, int length2, int head2, bool isInt = false);

//private:
    double InterpLine(double xi,vector<double>&x,vector<double>&y);
};

class BeiDouSV:public SV{
public:
    //signal number:2,3,4;  8,9,10;  14,15,16;
    SignalData B1_2I,B1_2Q,B1_2X,B3_6I,B3_6Q,B3_6X,B2_7I,B2_7Q,B2_7X;
public:
    BeiDouSV(int id);
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
    GpsSV(int id);
    ~GpsSV();
    int DecodeSubFrame(uint32_t* dwrds);
    SignalData* SignalTable(int index);
};

class SvSys{
public:
    SysType type;
    vector<SV*> table,used;
    SvSys(SysType _type):type(_type){};
    void OpenClose(bool state){
        for(SV* sv:table)sv->open=state;
    }
};

class SvAll{
public:
    GNSS* gnss;
//    SvSys* collect[Nsys];
//    vector<int> sysIds;
    vector<SvSys*> sysAll;
    vector<SV*> svUsedAll;

//    GpsSV svGpss[Ngps];
//    BeiDouSV svBeiDous[Nbds];
public:
    SvAll();
    void SetOpen(bool bds,bool gps);
//    SvAll(GNSS* gnss);
    ~SvAll();
    void InitAlloc();
    void UpdateEphemeris(char * subFrame);
    SV* GetSv(SysType type, int ind){
        SvSys *sys =GetSys(type);
        if(ind>sys->table.size()+1){
            printf("svSys Data not found!\n");
            return nullptr;
        }
        return sys->table[ind];
    }
    SvSys* GetSys(SysType type){
        for(SvSys* sys:sysAll){
            if(sys->type==type)return sys;
        }
        printf("svSys Data not found!\n");
        return nullptr;
    }
    int AddUsed(SV *sv);
private:
};

#endif