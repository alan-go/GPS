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
    GnssTime timeSolver;
    Vector3d xyz, vxyz;
    Solution soltion;
    vector<SV*> svsForCalcu[Nsys];

    double tu[Nsys]={0};
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

    int nSat,nSys;

    int PrepareSVsData(vector<SV*> &_svsIn);
    int SelectSvsFromVisible(vector<SV*> &all);
    int ProcessRtkData();
    int UpdateSvsPosition(vector<SV*> &svs, GnssTime rt, int ephType);
    int SolvePosition(vector<SV*>svsForCalcu);
    int SolvePositionBeiDouGPS(vector<SV*>svsForCalcu);
    int SolvePositionCalman();

};

class PosSolverKalman:private PosSolver{
public:
    VectorXd cyclS,Pii;
    Matrix<double,6,6> Pxv;

    int PositionKalman(vector<SV*> _svsIn);
};

#endif