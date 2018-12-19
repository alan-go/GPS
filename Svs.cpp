#include "Svs.h"
#include "NtripRTK.h"
#include "EphemSp3.h"
#include "EphemBst.h"

SV::SV(int id):svId(id),I(0),T(0),isBeiDouGEO(false),elevationAngle(0),tsdt(0),trackCount(0){
    kal = Kalman(30,29,0);
    ephemSp3 = new EphemSp3(this);
    ephemBst = new EphemBst(this);
}
SV::SV(){}
SV::~SV(){}
BeiDouSV::BeiDouSV(int id):SV(id){type = SYS_BDS;}
BeiDouSV::~BeiDouSV(){}
GpsSV::GpsSV(int id):SV(id){type = SYS_GPS;}
GpsSV::~GpsSV(){}


void SvAll::SetOpen(bool bds, bool gps) {
    GetSys(SYS_BDS)->OpenClose(bds);
    GetSys(SYS_GPS)->OpenClose(gps);
}
SvAll::SvAll(){}

void SvAll::InitAlloc() {
    SvSys* sBds = new SvSys(SYS_BDS);
    SvSys* sGps = new SvSys(SYS_GPS);
    for (int i = 0; i < Nbds; ++i) {
        SV* sv = new BeiDouSV(i+1);
        sBds->table.push_back(sv);
        sv->isBeiDouGEO = i<5?true:false;
        sv->ephemBst->ion = &ionKlob;
        svAll.push_back(sv);
    }
    sysAll.push_back(sBds);

    for(int i=0;i<Ngps;i++){
        SV* sv = new GpsSV(i+1);
        sGps->table.push_back(sv);
        sv->ephemBst->ion = &ionKlob;
        svAll.push_back(sv);
    }
    sysAll.push_back(sGps);
}

//todo get usedSys
int SvAll::AddToUsed(SV *sv) {
    svUsed.push_back(sv);
    GetUsedSys(sv->type)->table.push_back(sv);
}
//
int SvAll::RemoveSvFromUsed(SV *sv){
    vector<SV*> temp;
    for(SV* svi:svUsed){
        if(svi==sv)continue;
        temp.push_back(svi);
    }
    svUsed.swap(temp);
    UpdateSVs();
}
SvAll::~SvAll(){}


bool SV::CheckEphemStates(GnssTime time,int ephemType) {
    ephemOkBst = ephemBst->Available(time);
    ephemOkSp3 = ephemSp3->Available(time);
    switch (ephemType){
        case 0:
            return ephemOkBst;
        case 1:
            return ephemOkSp3;
        default:
            return 0;
    }
}



bool SV::AddMmeasure(Measure *mesr) {
    if(!measureDat.empty()){
        Measure* latest = measureDat.front();
        mesr->cycle = latest->cycle;
        mesr->cycleP = latest->cycleP;
        if(mesr->time-latest->time>1.5)trackCount=0;
        else trackCount++;
    }
    measureDat.push_front(mesr);
};

void SV::PrintInfo(int printType) {
    //1:position and tsDelta
    //2:bst ephemeris
    switch (printType){
        case 0:
            printf("sv:%d,%02d  %s\n", type, svId,tip);
            break;
        case 1:
//            cout<<"+++++++SvPosition:"<<type<<","<<svId<<endl;
//            cout<<position<<endl;
            printf("+++SvPosition:%d,%d === %10f,%10f,%10f\n",type,svId,xyz(0),xyz(1),xyz(2));
            cout<<"norm="<<xyz.norm()<<endl;
            printf("+++svLLA === %lf, %lf, %lf\n",lla(0)*180/GPS_PI,lla(1)*180/GPS_PI,lla(2));
//            printf("+++TGD === %.3f\n",TGD1*1e9);
//
//            cout<<"tsDelta === "<<tsdt*Light_speed<<" ,a0"<<a0<<endl;
            break;
        case 2:
            printf("\n%d,%d\n",type,svId);
//            printf("%d,\t%.10e,\t%.10e,\t%.10e,\t\n",SOW,a0,a1,a2);
//            printf("%d,\t%.10e,\t%.10e,\t%.10e,\t\n",orbit.IODE,orbit.Crs,orbit.dtn,orbit.M0);
//            printf("%.10e,\t%.10e,\t%.10e,\t%.10e,\t\n",orbit.Cuc,orbit.e,orbit.Cus,orbit.sqrtA);
//            printf("%d,\t%.10e,\t%.10e,\t%.10e,\t\n",orbit.toe,orbit.Cic,orbit.Omega0,orbit.Cis);
//            printf("%.10e,\t%.10e,\t%.10e,\t%.10e,\t\n",orbit.i0,orbit.Crc,orbit.omega,orbit.OmegaDot);
//            printf("%.10e,\t%d,\t%d,\t%d,\t\n",orbit.IDOT,0,WN,0);
//            printf("%d,\t%d,\t%.10e,\t%d,\t\n",0,SatH1,TGD1,IODC);
            break;
        default:
            break;
    };
}

void SV::FPrintInfo(int printType){
    if(NULL==fpLog){
    char name[32];
    sprintf(name,"../log/SV/%d_%02d.txt",type,svId);
    fpLog = fopen(name,"w");
    }
    Measure *ms = measureDat.front();
    Measure *ms1 = measureDat[1];
    double dt = ms->time.tow-(ms1->time.tow);
    switch(printType){
        case 0:

            break;
        case 1:
            fprintf(fpLog,"%.4f,\t%.4f,%.4f,\t",ms->time.tow,ms->prMes,ms->stdevPr);//1,2
            fprintf(fpLog,"%.4f,%.4f,\t",ms->cpMes,ms->stdevCp);//3,4
            fprintf(fpLog,"%.4f,%.4f,\t",ms->doMes,ms->stdevDo);//5,6
            fprintf(fpLog,"%.4f,%d,\t",ms->lockTime/1000,ms->trkStat);//7,8
            fprintf(fpLog,"%.4f,%.4f,%.4f\n",ms->prMes-ms->cpMes,(ms->prMes-ms1->prMes)/dt,(ms->cpMes-ms1->cpMes)/dt);//9,10,11
            break;
        case 2:
            break;
    }
}

void SvAll::UpdateEphemeris(char *subFrame) {
    char* playload = subFrame + 6;
    uint8_t gnssId = *(uint8_t*)(playload);
    uint8_t svId = *(uint8_t*)(playload + 1);
    uint8_t numWords = *(uint8_t*)(playload+4);
    printf("Update subframe;;gnssid:%d,svid:%d\n",gnssId,svId);
    char* tmp = playload+8;
    if(10==numWords){
        uint32_t dwrds[10];
        for(int i=0;i<10;i++) {
            dwrds[i] = *(uint32_t*)(tmp+4*i);
        }
        SV* sv = GetSv(SysType(gnssId),svId);
        if(nullptr!=sv)sv->ephemBst->DecodeSubFrame(dwrds);
    }
}

int SV::CalcuelEvationAzimuth(Vector3d pos, Vector3d poslla) {
    Vector3d dtenu,dtxyz = xyzR - pos;
    XYZ2ENU(dtxyz,poslla,dtenu);
    elevationAngle = asin(dtenu(2)/dtenu.norm());
    azimuthAngle = atan2(dtenu(0),dtenu(1));
    if(azimuthAngle<0)azimuthAngle+=2*GPS_PI;
}

int SV::CalcuInoshphere(double elev, double azim,Vector3d LLA,double time) {
    Ionosphere *ion = ephemBst->ion;
    double temp0 = Earth_a/(Earth_a+375000);
    double temp1 = temp0*cos(elev);
//    printf("ai,bi %e,%ef,%ef,%ef,\t %e,%e,%e,%e,\n",ion->a0,ino.a1,ino.a2,ino.a3,ino.b0,ino.b1,ino.b2,ino.b3 );
    double psi = GPS_PI/2 - elev - asin(temp1);
    double phyM = asin(sin(LLA(0))*cos(psi) + cos(LLA(0))*sin(psi)*cos(azim));
    double lambdaM = LLA(1)+asin(sin(psi)*sin(azim)/cos(phyM));
    double t = time + lambdaM*43200/GPS_PI;
    t = fmod(t,86400);
    if(t<0)t+=86400;

    double phyMpi = phyM/GPS_PI;
    double A2 = ion->a0 + phyMpi * (ion->a1 + phyMpi * (ion->a2 + phyMpi * ion->a3));
    if(A2<0)A2 = 0;
    double A4 = ion->b0 + phyMpi * (ion->b1 + phyMpi * (ion->b2 + phyMpi * ion->b3));
    if(A4<72000)A4 = 72000;
    if(A4>172800)A4 = 172800;

    double Iz_ = 5e-9;
    if(abs(t - 50400) < A4/4)
        Iz_ += A2*cos(2*GPS_PI*(t-50400)/A4);

    I = Iz_ / sqrt(1 - temp1*temp1);
    I *= Light_speed;
    return 0;
}

int SV::CalcuTroposhphere(double elev, double azim) {
    if(elev<0.1116)
        T=20;
    else
        T=2.47/(sin(elev)+0.0121);
    return 0;
}

int SV::CorrectIT(Vector3d xyz,Vector3d lla,double time) {
    CalcuelEvationAzimuth(xyz,lla);
//    printf("%d,%d..elevation = %lf, azim = %lf\n",type,svId,elevationAngle,azimuthAngle);
    CalcuTroposhphere(elevationAngle,azimuthAngle);
    CalcuInoshphere(elevationAngle,azimuthAngle,lla,time);
}



bool SV::ElevGood() {
    if(0==elevationAngle)
        return 1;
    elevGood = elevationAngle>0.17;
    if(!elevGood){
        sprintf(tip,"elevTooLow:%.1f",elevationAngle);
    }
    return elevGood;
}

bool SV::MeasureGood() {
    Measure *ms0 = measureDat.front();
    measureGood = ms0->prMes<45e6&&ms0->prMes>15e6;

    if(!measureGood){
        sprintf(tip,"prError:%.1f",ms0->prMes);
    }
    if (trackCount<5) {
        sprintf(tip,"trackCount %d<5", trackCount);
        return 0;
    }
//    if(ms0->trkStat<3)return 0;

    return measureGood;
}

double SV::InterpRtkData(double time, int sigInd) {
    int InterpLength = 5;

    auto ite = rtkData.begin();
    double dt = 100;

    for (; ite < rtkData.end() - 5; ite++) {
        dt = time - (*ite)->rtktime;
        if(dt>0)break;
    }
    if(dt>2)
    {
        sprintf(tip," %drtk,cloest=%.1f",rtkData.size(),dt);
        return 0;
    }
//    printf("time: %lf,%lf\n", time,(*ite)->rtktime);

    vector<double> times, prMes, cpMes;
    for(int i=0;i<InterpLength;i++){
        MSM4data* data = *(ite+i);
        MSM4Cell* cell = &data->sigData[sigInd];
        if(cell->cpMes!=0){
            times.push_back(data->rtktime);
            prMes.push_back(cell->prMes);
            cpMes.push_back(cell->cpMes);
//            printf("t=%f,pr=%f,cp = %f\n",data->rtktime,cell->prMes,cell->cpMes);
        }

    }
    if(times.size()<InterpLength-2)
    {
        printf("interp time line size = %d, not ok\n",times.size());
        return 0;
    }
    prInterp[sigInd] = InterpLine(time,times,prMes);
    cpInterp[sigInd] = InterpLine(time,times,cpMes);
//    printf("Interp:t=%f,pr=%f,cp = %f\n",time,prInterp[sigInd],cpInterp[sigInd]);

    return prInterp[sigInd];
}

double SV::InterpLine(double xi, vector<double> &x, vector<double> &y) {
    if(x.size()!=y.size()) return 0;
    int N = x.size();
    MatrixXd x2(N,2);
    MatrixXd y1(N,1);
    for(int i=0;i<N;i++){
        x2(i,0) = 1;
        x2(i,1) = x[i]-x[0];
        y1(i) = y[i]-y[0];
//        printf("xi,yi = %f,%f\n",x2(i,1),y1(i));
    }
    MatrixXd x2T = x2.transpose();
    Vector2d ab = (x2T*x2).inverse()*(x2T*y1);
//    printf("ab = %f,%f\n",ab(0),ab(1));
    return ab(0) + ab(1)*(xi-x[0]) + y[0];
}

bool SV::IsMaskOn() {
    switch (type){
        case SYS_BDS:
            break;
    }
}
bool SV::RotateECEF(double tt) {
    EarthRotate(xyz,xyzR,tt);
    EarthRotate(vxyz,vxyzR,tt);
    XYZ2LLA(xyz,lla);
    XYZ2LLA(xyzR,llaR);
}

double SV::InterpMeasere(int len, int power, int begin) {
    if(measureDat.size()<(len+begin))
        return -1;
    MatrixXd G(len,power+1),W(len,len);
    VectorXd Ai(power+1),b(len);
    W.fill(0);
    double t0 = measureDat[begin]->time.tow;
    double dt1 = t0-measureDat[begin+1]->time.tow;
    if(dt1>3){
//        measureDat[0]->cycle=0;
        return -1;
    }
    for(int i=begin;i<begin+len;i++){
        Measure* ms = measureDat[i];
        double dt = t0-ms->time.tow;
        for(int n=0;n<=power;n++)G(i,n)=pow(dt,n);
        W(i,i)=1/pow(ms->stdevDo,2);
        b(i)=ms->doMes;
    }
    MatrixXd WG = W*G;
//    MatrixXd WGT = WG.transpose();
    Ai = WG.colPivHouseholderQr().solve(W*b);
    //积分
    double dDopler =0;
    for(int n=0;n<=power;n++){
        dDopler+=Ai(n)*pow(dt1,n+1)/(n+1);
    }


    /////////////debug
    double lambda = GetFreq(type,1,1);
    Measure* ms0 = measureDat[begin];
    Measure* ms1 = measureDat[begin+1];
    double dcp = ms0->cpMes-ms1->cpMes;
    double dpr = ms0->prMes-ms1->prMes;
//    dDopler = (ms0->doMes+ms1->doMes)*dt1/2;
//    if(abs(dcp-dDopler)>500)  {
//        return -1;
//    }
    double Cycle = ((dcp-dDopler)/lambda);
    double dCycle = fmod(Cycle+1e6,1);
    if(dCycle>0.5)dCycle-=1;
    double temp = (Cycle-dCycle);
    ms0->cycle+=temp;
    double pr_cpC0 = ms0->prMes-ms0->cpMes;
    double pr_cpC_ = ms0->prMes-ms0->cpMes-(ms0->cycle)*lambda;
    double pr_cpC1 = ms0->prMes-ms0->cpMes+(ms0->cycle)*lambda;
    fprintf(fpLog,"%.4f,\t",ms0->time.tow);//0,
    fprintf(fpLog,"%f,%f,%f,%f,%f,    ",dpr,dcp,dDopler,Cycle,dCycle);//12345
    fprintf(fpLog,"%f,%f,%f,%f,   ",(dpr)/dt1,(dDopler)/dt1,ms0->stdevPr,ms0->stdevDo);//6789
    fprintf(fpLog,"%f,%f,%f\n",pr_cpC0,pr_cpC_,pr_cpC1);//10 11
}

double SV::SmoothPr(int len, int begin) {
    double lambda = GetFreq(type,1,1);
    Measure* ms0 = measureDat[begin];
    Measure* ms1 = measureDat[begin+1];
    double t0 = ms0->time.tow,t1 = ms1->time.tow;
    double weight0 = 1/pow(ms0->stdevDo,2),weight1 = 1/pow(ms1->stdevDo,2);
    double dt01 = t0-t1;
    if(dt01<1){ms0->trackTime=ms1->trackTime+dt01;}
    else { return -1;}

    double dDoppler = (weight0*ms0->doMes+weight1*ms1->doMes)/(weight0+weight1)*dt01;
    double dcp = ms0->cpMes-ms1->cpMes;
    double dDoCp = dcp-dDoppler;
    ms0->cycleSlip = round(dDoCp/lambda);
    ms0->cycleRes = dDoCp/lambda-ms0->cycleSlip;

    int N = 0;
    double cycleTemp=0,diffSum=0,diffMean=0;
    for(int i=begin;;i++,N++){
        Measure* ms = measureDat[i];
        if(0==ms->trackTime||ms->stdevDo>0.5)break;
        if(t0-ms->time.tow>len)break;
        cycleTemp+=ms->cycleSlip;
        double diff = ms->cpMes-ms->prMes-cycleTemp*lambda;
        diffMean = (diff+diffMean*N)/(N+1);
    }

    double stdSum=0;
    for(int i=begin;i<N;i++){
        Measure* ms = measureDat[i];
        cycleTemp+=ms->cycleSlip;
        double diff = ms->cpMes-ms->prMes-cycleTemp*lambda;
        double std = pow(diff-diffMean,2);
        stdSum+=std;
    }

}


double SV::SmoothKalman(int len, int begin) {
    Measure* ms0 = measureDat[begin];
    Measure* ms1 = measureDat[begin+1];
    double t0 = ms0->time.tow,t1 = ms1->time.tow;
    double dt01 = t0-t1;

    if(dt01>1.5){
        kal.state=0;
        return -1;
    }
    double weight0 = 1/pow(ms0->stdevDo,2),weight1 = 1/pow(ms1->stdevDo,2);
    if(dt01<1){ms0->trackTime=ms1->trackTime+dt01;}
    else { return -1;}
    ms0->dDoppler = (weight0*ms0->doMes+weight1*ms1->doMes)/(weight0+weight1)*dt01;

    kal.state++;
    int N = kal.N,M=kal.M;
    if(kal.state<kal.N)return -2;
    if(kal.state==kal.N){
        kal.Init();
        for(int ni=0;ni<N;ni++)kal.x(ni)=measureDat[ni]->prMes;
    }

    for(int ni=0;ni<N;ni++){
        Measure* msi = measureDat[ni];
        kal.y(ni)=msi->prMes;
        kal.Rmm(ni,ni) = pow(msi->stdevPr,2);
        if(N-1==ni)continue;
        kal.y(ni+N)=msi->dDoppler;
        kal.Rmm(ni+N,ni+N)=pow(msi->stdevDo,2);
    }

    kal.x.tail(N-1)=kal.x.head(N-1);
    kal.x(0)+=ms0->dDoppler;
    kal.Hmn.fill(0);
    kal.Hmn.block(0,0,N,N)=MatrixXd::Identity(N,N);
    kal.Hmn.block(N,0,N-1,N-1)+=MatrixXd::Identity(N-1,N-1);
    kal.Hmn.block(N,1,N-1,N-1)-=MatrixXd::Identity(N-1,N-1);
//    cout<<"Rn:"<<endl<<kal.Rmm<<endl;
//    cout<<"P:"<<endl<<kal.Pnn<<endl;

    fprintf(fpLog,"%.4f,\t %.10f",ms0->time.tow,ms0->prMes);//0,
//    fprintf(fpLog,"%.4f,\t %.10f",ms0->time.tow,kal.x(0));//0,

    kal.Rectify();
    fprintf(fpLog,",%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",kal.x(0),(kal.x(0)-kal.x(1))/dt01,ms0->dDoppler/dt01);//0,
}

double SV::SmoothKalman0() {
    Measure* ms0 = measureDat[0];
    ms0->prCor = ms0->prMes;
    ms0->stdPrCor = ms0->stdevPr;
    if(measureDat.size()<10)
        return -1;
    Measure* ms1 = measureDat[1];
    double t0 = ms0->time.tow,t1 = ms1->time.tow;
    double dt01 = t0-t1;

    if(dt01>1.5){
        kal.state=0;
        return -1;
    }
    double weight0 = 1/pow(ms0->stdevDo,2),weight1 = 1/pow(ms1->stdevDo,2);
    if(dt01<1){ms0->trackTime=ms1->trackTime+dt01;}
    else { return -1;}
    ms0->dDoppler = (weight0*ms0->doMes+weight1*ms1->doMes)/(weight0+weight1)*dt01;
    double dcp = ms0->cpMes-ms1->cpMes;

    double lambda = GetFreq(type,1,1);
    double dDoCp = dcp-ms0->dDoppler;
    ms0->cycleSlip = round(dDoCp/lambda);
    ms0->cycleRes = dDoCp/lambda-ms0->cycleSlip;

    kal.state++;
    int N = kal.N,M=kal.M;
    if(kal.state<kal.N)return -2;

    kal.Pnn.fill(0);
    for(int ni=0;ni<N;ni++){
        Measure* msi = measureDat[ni];
        kal.x(ni)=msi->prMes;
        kal.Pnn(ni,ni)=pow(msi->stdevPr,2);

        if(N-1==ni)continue;
        kal.y(ni)=msi->dDoppler;
        kal.Rmm(ni,ni)=pow(msi->stdevDo,2);
    }


    kal.Hmn.fill(0);
    kal.Hmn.block(0,0,N-1,N-1)+=MatrixXd::Identity(N-1,N-1);
    kal.Hmn.block(0,1,N-1,N-1)-=MatrixXd::Identity(N-1,N-1);
//    cout<<"Rn:"<<endl<<kal.Rmm<<endl;
//    cout<<"P:"<<endl<<kal.Pnn<<endl;

    fprintf(fpLog,"%.4f,\t %.10f",ms0->time.tow,ms0->prMes);//0,

    kal.Rectify();
    ms0->prCor = kal.x(0);
    ms0->stdPrCor = sqrt(kal.Pnn(0,0));
//    cout<<"Pafter:"<<endl<<kal.Pnn<<endl;
    fprintf(fpLog,",%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",kal.x(0),ms0->stdPrCor,ms0->prMes-kal.x(0));//0,
}
double SV::DetectCycleSlip() {
    if(measureDat.size()<2)
        return -1;
    Measure* ms0 = measureDat[0];
    Measure* ms1 = measureDat[1];
    double t0 = ms0->time.tow,t1 = ms1->time.tow;
    double dt01 = t0-t1;

    if(dt01>1.2){
        ms0->trackTime=0;
        return -1;
    } else{
        ms0->trackTime=ms1->trackTime+dt01;
    }
    double weight0 = 1/pow(ms0->stdevDo,2),weight1 = 1/pow(ms1->stdevDo,2);
    ms0->dDoppler = (weight0*ms0->doMes+weight1*ms1->doMes)/(weight0+weight1)*dt01;
    double dcp = ms0->cpMes-ms1->cpMes;

    double lambda = GetFreq(type,1,1);
    double dDoCp = dcp-ms0->dDoppler;
    ms0->cycleSlip = round(dDoCp/lambda);
    ms0->cycleRes = dDoCp/lambda-ms0->cycleSlip;
    //调参？
//    ms0->cycleSlipQ = pow(10*dt01*dt01/lambda,2);
    ms0->cycleSlipQ = pow(50/lambda*ms0->stdevDo*dt01*dt01,2);
}
