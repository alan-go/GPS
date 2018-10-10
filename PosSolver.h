#ifndef POSSOLVER_H
#define POSSOLVER_H

#include "Svs.h"
#include "NtripRTK.h"
#include "CommonInclude.h"


class PosSolver{
public:
    SvAll svsBox;
    GNSS* gnss;
    NtripRTK *rtk;
    double PDOP;
//    double rcvtow;
    GnssTime timeSolver;
    Vector3d xyz, vxyz;
    Solution soltion;
    vector<SV*> svsForCalcu[Nsys];
//    double tuBds, tuGps;
    double tu[Nsys]={0},tui;
    int numMeas;
public:
    PosSolver();
    PosSolver(SvAll svs, NtripRTK *rtk, GNSS *gnss);
    ~PosSolver();
    int PositionSingle(vector<SV*> _svsIn);
    int PositionRtk();
    int PositionRtk(vector<SV*> _svsIn);
    int PositionRtk2();
    int PositionRtkKalman();
    int MakeGGA(char *gga,Vector3d lla,GnssTime gpsTime);

private:
    int sysEachNum[Nsys],sysCount;
//    int N,N_n;
    int nSat,nSys;
private:
    int PrepareSVsData(vector<SV*> &_svsIn);
    int SelectSvsFromVisible(vector<SV*> &all);
    int ProcessRtkData();
    int UpdateSvsPosition(vector<SV*> &svs, GnssTime rt, int ephType);
    int SolvePosition(vector<SV*>svsForCalcu);
    int SolvePositionBeiDouGPS(vector<SV*>svsForCalcu);
    int SolvePositionCalman();

};

#endif