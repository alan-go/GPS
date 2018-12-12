//
// Created by root on 8/31/18.
//
#include "GNSS.h"


#ifndef GPS_EPHEMSP3_H
#define GPS_EPHEMSP3_H



struct Sp3Cell{
    GnssTime time;
    Eigen::Vector3d pxyz, vxyz;
    double ts, tsDrift;
    Sp3Cell():ts(0),tsDrift(0){
        pxyz<<0,0,0;
        vxyz<<0,0,0;
    }
};
struct TsCell{
    GnssTime time;
    double ts{0};
    TsCell(GnssTime time_,double ts_):time(time_),ts(ts_){}
};

class EphemSp3 {
public:
    SV* sv{nullptr};
    vector<Sp3Cell> xyzList;
    vector<TsCell> tsList;
    GnssTime timeHead, timeEnd;
    double dt = 15;
public:
    EphemSp3();
    EphemSp3(SV* sv);
    static int ReadSp3File(string fileName, SvAll &svs);
    static int ReadClkFile(string fileName, SvAll &svs);
    static int ReadSp3s(string name, SvAll &svs);
    static void* GetSp3Thread(void* _gnss);
    static SysType code2sys(char code);
    int CalcuECEF(GnssTime interpTime);
    int CalcuTs(GnssTime ts0);
    int interpECEF(GnssTime interpTime,Vector3d &xyz);
    bool Available(GnssTime time);
};




#endif //GPS_EPHEMSP3_H
