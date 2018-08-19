#ifndef POSSOLVER_H
#define POSSOLVER_H

#include "SVs.h"
#include "NtripRTK.h"
#include <Eigen/Core>
#include <Eigen/Dense>

class PosSolver{
public:
    SVs svs;
    GNSS* gnss;
    char raw[1024];
    NtripRTK *rtk;
    double rcvtow;
    Vector3d xyz, LLA;
    double tu, tuBeiDou, tuGps;
    int numMeas;
public:
    PosSolver(SVs svs, NtripRTK *rtk, GNSS *gnss);
    ~PosSolver();
    int PositionSingle();
    int PositionRtk();
    static int XYZ2LLA(Vector3d XYZ,Vector3d &LLA);

private:
    vector<SV*>visibleSvs;
private:
    int ReadVisibalSvsRaw(SVs svs, vector<SV*> &svVisable, char *raw);
    int SelectSvsFromVisible(vector<SV*> &all,vector<SV*> &select,int &bds,int &gps);
    int SolvePosition(vector<SV*>svsForCalcu);
    int SolvePositionBeiDouGPS(vector<SV*>svsForCalcu);
    int SolvePositionCalman();
    int InosphereCorrect();
    int TroposphereCorrect();
};

#endif