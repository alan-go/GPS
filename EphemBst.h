//
// Created by alan on 18-12-12.
//

#ifndef GPS_EPHEMBST_H
#define GPS_EPHEMBST_H

#include "GNSS.h"



class EphemBst {
public:
struct Orbit{
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
    Orbit():toe(0),sqrtA(0),e(0),i0(0),Omega0(0),M0(0), Cus(0),Cuc(0),Cis(0),
            Cic(0),Crs(0),Crc(0),dtn(0),omega(0),IDOT(0),OmegaDot(0){}
};

struct Clock{
    double toc{0},a0{0},a1{0},a2{0},a3{0};
    double TGD1{0},TGD2{0};
    uint32_t a1High{0},a1Low{0};
    //week number; second of week;
    uint32_t SOW{0},WN{0};
    //时钟数据龄期(AODC/IODC),星历数据龄期(AODE/IODE)
    uint32_t AODC,AODE;
};
public:
    SV* sv{nullptr};
    Orbit orbit;
    Ionosphere *ion;
    Clock clk;

    int8_t framOk[10];
    //health,用户距离精度指数
    uint32_t SatH1{1},URAI{15};
public:
    EphemBst();
    EphemBst(SV* sv_);
    static int ReadFromFile(string fileName, SvAll &svs);
    static void* DownLoadThread(void* _gnss);
    int CalcuECEF(GnssTime time);
    int CalcuTs(GnssTime ts0);
    bool Available(GnssTime time);

    int DecodeSubFrame(uint32_t* dwrds);

private:

    int DecodeBdsD1(uint32_t* dwrds);
    int DecodeBdsD2Frame1(uint32_t* dwrds);
    int DecodeGps(uint32_t* dwrds);
    //head 指32bit()中的头bit（范围：1-32）
    inline uint32_t Read1Word(uint32_t word, int length, int head, bool isInt = false);
    inline uint32_t Read2Word(uint32_t word0,int length0, int head0,
                              uint32_t word1, int length1, int head1, bool isInt = false);
    inline uint32_t Read3Word(uint32_t word0,int length0, int head0,
                              uint32_t word1, int length1, int head1, uint32_t word2, int length2, int head2, bool isInt = false);

};


#endif //GPS_EPHEMBST_H
