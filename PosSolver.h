#ifndef POSSOLVER_H
#define POSSOLVER_H

#include "SVs.h"
#include "NtripRTK.h"
#include "CommonInclude.h"

struct MeasData{
    GnssTime time;
    double prMes, cpMes, doMes;
};

class PosSolver{
public:
    SVs *svs;
    GNSS* gnss;
    char raw[1024];
    NtripRTK *rtk;
    double rcvtow;
    GnssTime rTime;
    Vector3d xyz, LLA, vxyz;
    double tu, tuBeiDou, tuGps;
    int numMeas;
public:
    PosSolver();
    PosSolver(SVs *svs, NtripRTK *rtk, GNSS *gnss);
    ~PosSolver();
    int PositionSingle();
    int PositionRtk();
    int PositionRtk2();
    int PositionRtkKalman();
    int MakeGGA(char *gga,Vector3d lla,GnssTime gpsTime);

private:
    vector<SV*>visibleSvs;
    int numOfSys[Nsys],nsysUsed;
private:
    int PrepareSVsData(vector<SV*> *svsOut);
    int ReadVisibalSvsRaw(SVs *svs, vector<SV*> &svVisable, char *raw);
    int SelectSvsFromVisible(vector<SV*> &all,vector<SV*> *select);
    int ProcessRtkData(vector<SV *> *select);
    int UpdateSvsPosition(vector<SV*> &svs, GnssTime rt, int ephType);
    int SolvePosition(vector<SV*>svsForCalcu);
    int SolvePositionBeiDouGPS(vector<SV*>svsForCalcu);
    int SolvePositionCalman();

};

#endif