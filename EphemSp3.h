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

class EphemSp3 {
public:
    vector<Sp3Cell> records;
    GnssTime timeHead, timeEnd;
    double dt = 15;
public:
    EphemSp3();
    static int ReadSp3File(string fileName, SVs *svs);
    static SV::SysType code2sys(char code);
    static Sp3Cell InterpECEF(vector<Sp3Cell> &list, GnssTime interpTime);
};


#endif //GPS_EPHEMSP3_H
