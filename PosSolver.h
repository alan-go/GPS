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
    int CalcuPosition();
    static int XYZ2LLA(Vector3d XYZ,Vector3d &LLA);

private:
    vector<SV*>visibleSvs, SvsForCalcu;
private:
    int PrepareData(SVs svs, char *raw);
    int SolvePosition();
    int SolvePositionBeiDouGPS();
    int SolvePositionCalman();
    int InosphereCorrect();
    int TroposphereCorrect();
};

#endif