//
// Created by root on 8/1/18.
//

#include "GNSS.h"
#include <stdio.h>
using namespace std;
using namespace Eigen;



int main(){
    printf("\nin testing\n");
    GNSS gnss;
    gnss.log = fopen("../log/log.txt","a+");

//    gnss.useGPS = false;
    gnss.useBeiDou = false;
//    gnss.useQianXun = false;


    SV sv[6];

    sv[0].position = Vector3d(10,20,8);
    sv[1].position = Vector3d(5,7,6);
    sv[2].position = Vector3d(12,5,3);
    sv[3].position = Vector3d(-2,4,0);
    sv[4].position = Vector3d(0,8,12);
    sv[5].position = Vector3d(3,-2,-3);

    sv[0].prMes = 120.73644;
    sv[1].prMes = 107.071068;
    sv[2].prMes = 111.40175;
    sv[3].prMes = 104.690416;
    sv[4].prMes = 110.86278;
    sv[5].prMes = 107.483315;

    vector<SV*> svs;
    SV* svtemp;
    for (int i = 0; i < 6; i++) {
        svtemp = &(sv[i]);
        svs.push_back(svtemp);
    }

    PosSolver slv;
    slv.xyz = Vector3d(0,0,0);
    slv.SolvePosition(svs);














    string serialPort = "/dev/ttyUSB0";
    gnss.StartGNSS(serialPort,115200);
    getchar();

    return 0;
}
