//
// Created by Haoliu on 8/1/18.
//

#include "GNSS.h"
#include <stdio.h>
#include "EphemSp3.h"
using namespace std;
using namespace Eigen;

void WriteSols(SolutionDeque sols,string saveName){
    if(sols.size()==0){
        printf("no data tobe write\n");
        return;
    }
    string path = "../log/"+saveName+".txt";
    FILE *fp = fopen(path.data(),"w");
    for(Solution sl:sols){
//       fprintf(fp,"%.4f,%.4f,%.4f,%.4f\n",sl.time.tow,sl.xyz(0),sl.xyz(1),sl.xyz(2));
        fprintf(fp,"%.4f,%.8f,%.8f,%.8f\n",sl.time.tod,sl.lla(0)*R2D,sl.lla(1)*R2D,sl.lla(2));
    }
    fclose(fp);
}


int main(int argc,char* argv[]){

    printf("\nin testing\n");
    GNSS *gnss = new GNSS();
    gnss->AddSerial(0,1,"/dev/ttyUSB0",115200,true);
    gnss->AddSerial(1,0,"/dev/ttyUSB1",115200,true);
    gnss->log = fopen("../log/log.txt","w");
    gnss->logDebug = fopen("../log/logDebug.txt","w");
    gnss->logTu = fopen("../log/logtu.txt","w");
    gnss->logSingle = fopen("../log/logSingle.txt","w");
    gnss->logSingleNew = fopen("../log/logSingleNew.txt","w");

    gnss->Init(0,0,1,1);

//    EphemSp3::ReadSp3s("/home/alan/projects/GPS/Sp3/hour20316_07",gnss->svsManager);
//    EphemSp3::ReadSp3s("/home/alan/projects/GPS/Sp3/hour20323_07",gnss->svsManager);
//    EphemSp3::ReadSp3s("/home/alan/projects/GPS/Sp3/hour20314_13",gnss->svsManager);
    FILE *fp;
    char name[128],dat[512],temp[256],tempc;
    string ss(argv[1]);
    ss = "1219_07_41";
    gnss->gpsWeek=2032;
    gnss->dow=3;

//    ss = "1215_07_56";
//    gnss->gpsWeek=2031;
//    gnss->dow=6;

    string ssData0 = "../data/device0_" + ss + ".data";
    string ssData1 = "../data/device1_" + ss + ".data";
    string ssRTK = "../data/" + ss + ".rtk";
    printf("open file name:%s\n",ssRTK.data());
    //ReadRTK data
    if(fp = fopen(ssRTK.data(),"rb")){
        int k = 0;
        while (!feof(fp)){

            dat[k] = fgetc(fp);
            if(0xd3==(u_char)dat[k-1]&&0==(dat[k]>>2)){
                gnss->rtkManager.ParaseRTK(dat,k-2);
                k=1;
                dat[0] = dat[k-1];
                dat[1] = dat[k];
                memset(dat+2,0, sizeof(dat)-2);
            }
            k++;
        }
        fclose(fp);
    } else printf("open rtk data failed \n");

    for(RefStation* rf:gnss->rtkManager.refs){
        printf("refStation %f,%d\n", rf->tow, rf->id);
        ShowV3(rf->pos,"pos");
    }

    //Reac raw1 measure data(ubx)
    if(fp = fopen(ssData1.data(),"rb")){
        while (128==fread(dat,1,128,fp)){
            gnss->GetSerial(1)->ScanSerialData(dat,128);
        }
        fclose(fp);
        WriteSols(gnss->GetSerial(1)->solRaw,ss+"xyzRAC");
    } else printf("open data failed \n");
    //Reac raw0 measure data(ubx)
    if(fp = fopen(ssData0.data(),"rb")){
        while (128==fread(dat,1,128,fp)){
            gnss->GetSerial(0)->ScanSerialData(dat,128);
        }
        fclose(fp);
        WriteSols(gnss->GetSerial(0)->solRaw,ss+"xyzUBX");
        WriteSols(gnss->solKalmans,ss+"xyzKAL");
        WriteSols(gnss->solKalDops,ss+"xyzKAL2");
        WriteSols(gnss->solRTKs,ss+"xyzRTK");
        WriteSols(gnss->solSingles,ss+"xyzSIG");

    } else printf("open data failed \n");
//    gnss->svsManager.GetSv(SYS_GPS,20)

//    getchar();
//system("python3 ../py/NewAnaPr.py");

    return 0;
}
