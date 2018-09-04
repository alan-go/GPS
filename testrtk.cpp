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
    gnss.log = fopen("../log/logqx.txt","w");

//    gnss.useGPS = false;
    gnss.useBeiDou = false;
//    gnss.useQianXun = false;
    gnss.Init(0,1,1,1);



    string serialPort = "/dev/ttyUSB0";
    gnss.StartGNSS(serialPort,115200);
    getchar();

    return 0;
}
