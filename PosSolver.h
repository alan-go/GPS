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
    double tu[Nsys],tuf;

    GnssTime timeSol,timeSolLast;//这个时间直接从接收机读出来的时间，有钟差
    double dts{0};
    Vector3d xyz, vxyz, lla;
    Solution soltion0,solSingle,solKalSigle,solKalDopp;
    vector<SV*> svsForCalcu[Nsys];

//    double tu[Nsys]={0};
    int nSat,nSys;

    int numMeas;
    Matrix<double,6,6> P66;
    Matrix<double,9,9> P99;
    Kalman kalSingle,kalRtk;
public:
    int InitKalman(GNSS *_gnss);
    int Init(GNSS *_gnss);
    PosSolver();
    PosSolver(SvAll svs, NtripRTK *rtk, GNSS *gnss);
    ~PosSolver();
    int PositionSingle(vector<SV*> _svsIn);
    int AnaMeasure();
    int PositionSingleNew(vector<SV*> _svsIn);
    int PositionRtk();
    int PositionRtk(vector<SV*> _svsIn);
    int PositionRtk2();
    int PositionRtkKalman();
    int MakeGGA(char *gga,Vector3d lla,GnssTime gpsTime);


    int PrepareSVsData(vector<SV*> &_svsIn);
    int SelectSvsFromVisible(vector<SV*> &all);
    int ProcessRtkData();
    int ProcessRtkData(GnssTime time,vector<SV*> svs);
    int UpdateSvsPosition(GnssTime rTime, int ephType);
    int SolvePosition(vector<SV*>svsForCalcu);
    int SolvePositionBeiDouGPS(vector<SV*>svsForCalcu);
    int SolvePositionCalman();
    int PositionKalman(vector<SV*> _svsIn);
    int PositionKalman2(vector<SV*> _svsIn);
    int AnaData(vector<SV*> _svsIn);

    int ResetKalSingle(int N,int M,vector<SV*> &svsIn,int L=0);
    int ResetKalRtk(int N,int M, Kalman & kal,int L=0);
    int PosKalSng(vector<SV *> _svsIn);

};



#endif