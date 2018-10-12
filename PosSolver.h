#ifndef POSSOLVER_H
#define POSSOLVER_H

#include "Svs.h"
#include "NtripRTK.h"
#include "CommonInclude.h"


class GNSS;
class PosSolver{
public:
    SvAll svsBox;
    GNSS* gnss;
    NtripRTK *rtk;
    double PDOP;

    GnssTime timeSolver;//这个时间直接从接收机读出来的时间，有钟差
    Vector3d xyz, vxyz, lla;
    Solution soltion;
    vector<SV*> svsForCalcu[Nsys];

    double tu[Nsys]={0};
    int numMeas;
    Matrix<double,6,6> Pxv;
public:
    int InitKalman(GNSS *_gnss);
    PosSolver();
    PosSolver(SvAll svs, NtripRTK *rtk, GNSS *gnss);
    ~PosSolver();
    int PositionSingle(vector<SV*> _svsIn);
    int PositionRtk();
    int PositionRtk(vector<SV*> _svsIn);
    int PositionRtk2();
    int PositionRtkKalman();
    int MakeGGA(char *gga,Vector3d lla,GnssTime gpsTime);

    int nSat,nSys;

    int PrepareSVsData(vector<SV*> &_svsIn);
    int SelectSvsFromVisible(vector<SV*> &all);
    int ProcessRtkData();
    int UpdateSvsPosition(vector<SV*> &svs, GnssTime rt, int ephType);
    int SolvePosition(vector<SV*>svsForCalcu);
    int SolvePositionBeiDouGPS(vector<SV*>svsForCalcu);
    int SolvePositionCalman();
    int PositionKalman(vector<SV*> _svsIn);

};



#endif