//
// Created by Haoliu on 8/1/18.
//

#include "GNSS.h"
#include <stdio.h>
using namespace std;
using namespace Eigen;

void WriteSols(SolutionDeque sols,string saveName){
    string path = "../log/"+saveName+".txt";
    FILE *fp = fopen(path.data(),"w");
    for(Solution sl:sols){
//       fprintf(fp,"%.4f,%.4f,%.4f,%.4f\n",sl.time.tow,sl.xyz(0),sl.xyz(1),sl.xyz(2));
        fprintf(fp,"%.4f,%.8f,%.8f,%.8f\n",sl.time.tod,sl.lla(0)*R2D,sl.lla(1)*R2D,sl.lla(2));
    }
    fclose(fp);
}


int main(){
    printf("\nin testing\n");
    GNSS *gnss = new GNSS();
    gnss->AddSerial(0,0,"/dev/ttyUSB0",115200,true,true);
    gnss->AddSerial(1,0,"/dev/ttyUSB1",115200,true,true);
    gnss->log = fopen("../log/log.txt","w");
    gnss->logDebug = fopen("../log/logDebug.txt","w");
    gnss->logTu = fopen("../log/logtu.txt","w");

    gnss->useQianXun = false;
    //ephem,qianxun,bds,gps
//    gnss->Init(0,0,0,1);
//    gnss->Init(0,0,1,0);
    gnss->Init(0,0,1,1);


    FILE *fp;
    char name[128],dat[512],temp[256],tempc;
//    string ss = "0921_13_02";
//    string ss = "1030_08_44";
    string ss = "1110_11_28";
//    string ss = "1110_13_07";
//        string ss = "0708-2.data";
//        string ss = "0823";
//        string ss = "0815-2";//soho novatal
//        scanf("%s",name);



    string ssData0 = "../data/device0_" + ss + ".data";
    string ssData1 = "../data/device1_" + ss + ".data";
    string ssRTK = "../data/" + ss + ".rtk";
    printf("open file name:%s\n",ssRTK.data());
    //ReadRTK data
    if(fp = fopen(ssRTK.data(),"rb")){
        int k = 0;
        while (!feof(fp)){
            dat[k++] = fgetc(fp);
            if(0xd3==(u_char)dat[k-2]&&0==(dat[k-1]>>2)){
                gnss->rtkManager.ParaseRTK(dat,k);
                k=2;
                dat[0] = dat[k-2];
                dat[1] = dat[k-1];
                memset(dat+2,0,510);
            }
        }
        fclose(fp);
    } else printf("open rtk data failed \n");

    //Reac raw1 measure data(ubx)
    if(fp = fopen(ssData1.data(),"rb")){
        while (128==fread(dat,1,128,fp)){
            gnss->GetSerial(1)->ScanSerialData(dat,128);
        }
        fclose(fp);
//        WriteSols(gnss->GetSerial(1)->solRaw,"xyzOf1030_08_44RAC");
    } else printf("open data failed \n");
    //Reac raw0 measure data(ubx)
    if(fp = fopen(ssData0.data(),"rb")){
        while (128==fread(dat,1,128,fp)){
            gnss->GetSerial(0)->ScanSerialData(dat,128);
        }
        fclose(fp);
//        WriteSols(gnss->GetSerial(0)->solRaw,"xyzOf1030_08_44UBX");
        WriteSols(gnss->solKalmans,"xyzOf1030_08_44KAL");
        WriteSols(gnss->solRTKs,"xyzOf1030_08_44RTK");
        WriteSols(gnss->solSingles,"xyzOf1030_08_44SIG");

    } else printf("open data failed \n");
//    gnss->svsManager.GetSv(SYS_GPS,20)

//    getchar();

    return 0;
}
