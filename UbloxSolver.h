#ifndef UBLOX_SOLVER_H
#define UBLOX_SOLVER_H

#include <iostream>
#include <cmath>

class SvInfo
{
    enum SvId{};
    struct Orbit{
        double sq_a, e,i0,omega0,w,M0;
        double Cus,Cuc,Cis,Cic,Crs,Crc,dtn,Omega,I;
        double toe,idoe;
    };
public:
    SvId id;
    Orbit orbit;
    double x,y,z;
    double prMesL1, snr, doppler;
public:
    bool CalcuECEF();
};

class UbloxSolver
{
public:
    double tow;
    double rx,ry,rz;//ECEF position of receiver
    double longitude,latitude,height;
    SvInfo GPSSVs[32];
    SvInfo BeiDouSVs[37];


public:
    UbloxSolver();
    ~UbloxSolver();
    void ParseRawData(const char* message);
    void ParseBstSubFrame(const char* message);
};

#endif //UBLOX_SOLVER_H
