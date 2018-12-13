//
// Created by root on 8/31/18.
//

#include "EphemSp3.h"

bool checkSp3ClkErp(string name0){
    bool result{true};
    if(0!=access((name0+".sp3").c_str(),0))result= false;
    if(0!=access((name0+".clk").c_str(),0))result= false;
    if(0!=access((name0+".erp").c_str(),0))result= false;
    return result;
}

EphemSp3::EphemSp3() {}
EphemSp3::EphemSp3(SV* sv_):sv(sv_) {}

//Download sp3 files from ftp://igs.gnsswhu.cn/pub/whu/MGEX/
void* EphemSp3::GetSp3Thread(void* _gnss){
    //whu一般整点后约半小时更新
    GNSS* gnss = (GNSS*)_gnss;
    string savePath = "../Sp3/";
    auto DownloadFtp = [](int week,string timeInd,string savePath)->bool{
        auto SysCmd = [](string saveName,string downPath){
            system(("curl -o "+saveName+downPath).c_str());
            if(0!=access(saveName.c_str(),0) ) return;
            system(("uncompress "+saveName).c_str());
        };
        string whuPath = " ftp://igs.gnsswhu.cn/pub/whu/MGEX/";
        char sweek[8];

        sprintf(sweek,"%d/",week);
        whuPath+=sweek;
        string nsp3 = timeInd+".sp3.Z";
        string nclk = timeInd+".clk.Z";
        string nerp = timeInd+".erp.Z";

        if(0!=access((savePath+timeInd+".sp3").c_str(),0))
            SysCmd(savePath+nsp3,whuPath+nsp3);
        if(0!=access((savePath+timeInd+".clk").c_str(),0))
            SysCmd(savePath+nclk,whuPath+nclk);
        if(0!=access((savePath+timeInd+".erp").c_str(),0))
            SysCmd(savePath+nerp,whuPath+nerp);

        if(0==access((savePath+"*.Z").c_str(),0))system(("rm "+savePath+"*.Z").c_str());
        return checkSp3ClkErp(savePath+timeInd);
    };
    mkdir(savePath.c_str(),S_IRWXU);
    struct tm *utcTime;
    time_t timeNow;
    while (!gnss->stop){
        timeNow = time(NULL);
        utcTime=gmtime(&timeNow);
        GnssTime time;
        time.epoch2time(utcTime);
        time.utc2gpst();
        int wday = floor(time.tow/86400);
        int hour = floor((time.tow-wday*86400)/3600);
        char timeInd[16];
        sprintf(timeInd,"hour%d%d_%02d",time.week,wday,hour);

        if(checkSp3ClkErp(savePath+timeInd)){
            printf("sp3 hour%d already have,sleep%d\n",hour, 60-utcTime->tm_sec);
            sleep(60-utcTime->tm_sec);
            continue;
        } else {
            if (DownloadFtp(time.week, string(timeInd), savePath)) {
                //Read in
                continue;
            } else {
                printf("sp3 hour%d not updated,trying hour%d,sleep%d\n",hour,hour-1, 60-utcTime->tm_sec);
                sprintf(timeInd, "hour%d%d_%02d", time.week, wday, hour - 1);
                if (checkSp3ClkErp(savePath + timeInd)) {
                    sleep(60 - utcTime->tm_sec);
                    continue;
                } else {
                    if (DownloadFtp(time.week, string(timeInd), savePath)) {
                        //read in data
                        continue;
                    } else{
                        printf("Warning!! Can not Download Sp3 from Internet,Tring againg...\n");
                        sleep(60);
                    }
                }
            }
        }
    }
}

int EphemSp3::ReadSp3File(string fileName, SvAll &svs) {
    FILE *fp;
    char buff[1024];
    char sp3Type;
    GnssTime sp3Time;
    int nsats;
    double neph,epht,timeAll;

    if (!(fp = fopen(fileName.data(), "r"))) {
        printf("sp3 file open failed. %s\n", fileName.data());
    }
    //Read Head
    for (int i=0;i<22;i++) {
        if (!fgets(buff,sizeof(buff),fp)) break;

        if (i==0) {
            sp3Type=buff[2];
            sp3Time = GnssTime(buff+3,28,0);
            neph = str2num(buff+32,7);
        } else if (1==i) {
            epht = str2num(buff+24,14);
            timeAll = epht * neph;
        } else if (2<=i&&i<=6) {
            if (i==2) {
                nsats=(int)str2num(buff+4,2);
            }
            for (int j=0,k=0;j<17&&k<nsats;j++) {
                SysType sys =code2sys(buff[9+3*j]);
                int prn=(int)str2num(buff+10+3*j,2);
            }
        }
    }
    //Read Body
    GnssTime ephTime;
    SV *sv;
    while (fgets(buff,sizeof(buff),fp)) {

        if (!strncmp(buff,"EOF",3)) break;

        if (buff[0]=='*'){
            ephTime = GnssTime(buff+3,28,0);
            if (ephTime.time == -1) {
                printf("sp3 invalid epoch %31.31s\n",buff);
                return -1;
            }
        } else{
            char dataType = buff[0];
            SysType tempSys =code2sys(buff[1]);
            int prn=(int)str2num(buff+2,2);
            sv = svs.GetSv(tempSys,prn-1);
            if(sv == nullptr)continue;

            Vector3d data(str2num(buff+4,14),str2num(buff+18,14),str2num(buff+32,14));
            Sp3Cell cell;
            cell.time = ephTime;
            if ('P'==dataType) {
                cell.pxyz = data*1000;
            }
            if ('V'==dataType) {
                cell.vxyz = data*1000;
            }
            sv->ephemSp3->xyzList.push_back(cell);
        }
    }
}

int EphemSp3::ReadClkFile(string fileName, SvAll &svs) {
    FILE *fp;
    char buff[128];
    if (!(fp = fopen(fileName.data(), "r"))) {
        printf("sp3 file open failed. %s\n", fileName.data());
    }
    while (fgets(buff, sizeof(buff), fp)) {
        if (!strncmp(buff, "EOF", 3)) break;
        if(strncmp(buff,"AS",2)) continue;
        SysType sys = code2sys(buff[3]);
        int prn=(int)str2num(buff+4,2);
        SV* sv = svs.GetSv(sys,prn-1);
        if(sv == nullptr)continue;
        //todo:Confirm: utc2gps??
        GnssTime time = GnssTime(buff+8,28,0);
        double ts = str2num(buff+40,20);
        sv->ephemSp3->tsList.push_back(TsCell(time,ts));
    }
}

int EphemSp3::ReadSp3s(string name, SvAll &svs) {
    //clear old data
    for(SV* sv:svs.svUsedAll){
        vector<TsCell> Ts0s;
        vector<Sp3Cell> xyz0s;
        Ts0s.swap(sv->ephemSp3->tsList);
        xyz0s.swap(sv->ephemSp3->xyzList);
    }
    //Read new data
    string nameSp3 = name+".sp3",nameClk = name+".clk";
    if(0==access(nameSp3.c_str(),0))ReadSp3File(nameSp3,svs);
    if(0==access(nameClk.c_str(),0))ReadClkFile(nameClk,svs);
}


SysType EphemSp3::code2sys(char code) {
    if (code=='G'||code==' ') return SYS_GPS;
    if (code=='C') return SYS_BDS; /* extension to sp3-c */
//    if (code=='R') return SYS_GLO;
//    if (code=='E') return SYS_GAL; /* extension to sp3-c */
//    if (code=='J') return SYS_QZS; /* extension to sp3-c */
//    if (code=='L') return SYS_LEO; /* extension to sp3-c */
    return SYS_NULL;
}

int EphemSp3::CalcuECEF(GnssTime interpTime) {
    Vector3d xyz0,xyz1;
    if(!Available(interpTime))return -1;
//    printf("interpTime= %d\n", interpTime.time);
    double dt = 1e-3;
    GnssTime time0 = interpTime,time1 = interpTime+dt;

    if(interpECEF(time0,xyz0)||interpECEF(time1,xyz1))return -1;

    //todo:satellite antenna offset correction

    sv->xyz = xyz0;
    sv->vxyz = (xyz1-xyz0)/dt;

    //todo:* relativistic effect correction */
    double dtRela = 2.0*sv->vxyz.dot(sv->vxyz)/Light_speed/Light_speed;
    sv->tsdt+=dtRela;
    //todo:     sv->tsDrift;
    return 0;
}

int EphemSp3::interpECEF(GnssTime interpTime,Vector3d &xyz) {
    double time[10],x[10],y[10],z[10];
//    int Nsp3 = xyzList.size();

    //calcu position
    int kt = FindTimeDate(xyzList,interpTime);
    GnssTime timeRef = xyzList[kt].time;
    for(int i=kt-4,j=0;i<kt+6;i++,j++){
        time[j]=xyzList[i].time-timeRef;
        x[j] = xyzList[i].pxyz(0);
        y[j] = xyzList[i].pxyz(1);
        z[j] = xyzList[i].pxyz(2);

    }
    double tt = interpTime-timeRef;
    xyz = Vector3d(lagrange(time,x,tt,10),lagrange(time,y,tt,10),lagrange(time,z,tt,10));

    return 0;
}

int EphemSp3::CalcuTs(GnssTime ts0) {
    double time[10],ts[10];
    if(!Available(ts0))return -1;

    int kt = FindTimeDate(tsList,ts0);
    GnssTime timeRef = tsList[kt].time;
    for(int i=kt-4,j=0;i<kt+6;i++,j++){
        time[j]=tsList[i].time-timeRef;
        ts[j] = tsList[i].ts;
    }

    sv->ts0=sv->ts= ts0;
    for(int i=0;i<3;i++){
        sv->tsdt = lagrange(time,ts,sv->ts-timeRef,10);
        sv->ts = ts0 - sv->tsdt;
    }
    return 0;
}

bool EphemSp3::Available(GnssTime time) {
    int Nts = tsList.size();
    if(tsList.size()<12 || time<tsList[5].time || tsList[Nts-6].time<time)return false;
    return true;
}