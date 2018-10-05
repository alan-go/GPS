#include "PosSolver.h"
#include "GNSS.h"
#include "EphemSp3.h"

PosSolver::PosSolver(){}

PosSolver::PosSolver(SVs *svs, NtripRTK *rtk, GNSS *gnss): svs(svs),rtk(rtk),gnss(gnss){
    memset(numOfSys,0,Nsys* sizeof(int));
    if(gnss->isPositioned){
        xyz = gnss->records.back().xyz;
        LLA = gnss->records.back().lla;
    } else{
        xyz = gnss->xyzDefault;
        LLA = gnss->llaDefault;
    }
    tu = gnss->tu;
    tuBeiDou = gnss->tuBeiDou;
    tuGps = gnss->tuGps;
}

PosSolver::~PosSolver(){}

int PosSolver::PrepareSVsData(vector<SV*> *svsOut) {
//    ReadVisibalSvsRaw(svs, visibleSvs, raw);
    ReadVisibalSvsRaw(gnss->svsManager, visibleSvs, raw);


////////////
//
//    if (visibleSvs[0]->measureRecord.size() > 20) {
//        for (int i = 0; i < visibleSvs.size(); i++) {
//            SV* sv = visibleSvs[i];
//            fprintf(gnss->log,"SV ID =  %d,%d\n", sv->type, sv->svId);
//            for (int j = 0; j < sv->measureRecord.size(); j++) {
//                SV::Measure* meas = &(sv->measureRecord[j]);
//                fprintf(gnss->log,"time,pr,cp =  %.5f, %.10f,%.10f,%.10f\n", meas->tow,meas->prMes,meas->cpMes,meas->doMes);
//            }
//            sv->measureRecord.clear();
//        }
//        fprintf(gnss->log,"\n\n\n\n\n\n\n\n\n\n");
//    }

//////////////

    SelectSvsFromVisible(visibleSvs,svsOut);

//////////////
//    for (int j = 0; j < svsOut.size(); ++j) {
//        SV *svTemp= svsOut[j];
//        if (svTemp->temp) {
//            svTemp->temp=0;
////            for (double t = rcvtow-360; t < rcvtow+360; t+=100) {
//            for (double t = 394352; t < 394960; t+=5) {
//                double timeForECEF = t;
//                if(SV::SYS_BDS == svTemp->type)timeForECEF-=14;
//                svTemp->CalcuTime(timeForECEF);
//                svTemp->CalcuECEF(svTemp->tsReal);
//                fprintf(gnss->log,"svsVis:time,%.4f,s%d%02d,pr,%.5f,prres,%.4f,norm-x-y-z,%.5f,%.5f,%.5f,%.5f\n",timeForECEF,
//                        svTemp->type,svTemp->svId,0.0,0.0,
//                        svTemp->position.norm(),svTemp->position(0),svTemp->position(1),svTemp->position(2));
//            }
//        }
//    }
//    return 0;
//////////////
    UpdateSvsPosition(svsOut[7], rTime, gnss->ephemType);

}

int PosSolver::PositionRtk2() {
    int sigInd = 1;
    //插值
    int result = 0;

    vector<SV*> svs0[Nsys], svs1;
    PrepareSVsData(svs0);
    printf("rcvtow = %lf\n",rcvtow);
    //暂时单模计算
    if(numOfSys[SYS_BDS]*numOfSys[SYS_GPS]){
        printf("NOT one sys mode ------------count 0\n");
//        return -1;
    }
    SV* svCectre;
    int svCentreInd = -2;
    double elevMax = -5;
//    gnss->rtkManager.mtxData.lock();
    printf(" positon rtk lock 0\n");
    printf("---------calcu N0 = %d\n",svs0[7].size());

    for(int i = 0,k=-1; i<svs0[7].size();i++){
        SV* sv = svs0[7][i];
        double time = rcvtow - tu/Light_speed;
        if(sv->type==SYS_BDS)time-=14;
        double pr = sv->InterpRtkData(time,sigInd);
        printf("--sv:%d,%d-pr time= %f, interp = %f\n",
               sv->type,sv->svId,time,pr);

        if(0.0!=pr){
            svs1.push_back(sv);
            k++;
            printf("angle = %.2f\n",sv->elevationAngle);
            if(sv->elevationAngle>elevMax){
                elevMax = sv->elevationAngle;
                svCectre = sv;
                svCentreInd = k;
            }
        }
    }
    printf("---------calcu N1 = %d, svCentreInd = %d\n",svs1.size(),svCentreInd);
    if(svs1.size()<4){
        printf("N1 not enough");
        return -1;
    }

    svs1.erase(svs1.begin()+svCentreInd);

    int N = svs1.size();
    printf("---------calcu N2 = %d\n",N);
    if(N<4){
//        gnss->rtkManager.mtxData.unlock();
        return -1;
    }
    MatrixXd pur(N,1);
    MatrixXd Gur(N,4);

    double pur0 = svCectre->prMes - svCectre->prInterp[sigInd];
    printf("pur0-base:%f - %f\n",svCectre->prMes ,svCectre->prInterp[sigInd]);
    Vector3d referBase = gnss->rtkManager.ECEF_XYZ;
    Vector3d referLLA;
    XYZ2LLA(referBase,referLLA);
//    printf("++++++++++++BaseLLA === %lf,%lf,%lf\n",referLLA(1)*180/GPS_PI,referLLA(0)*180/GPS_PI,referLLA(2));

    Vector3d Ir0 = svCectre->position - referBase;
    Ir0 = Ir0/Ir0.norm();

    printf("start 2222222222222SSSSSSSSSSSSSSSSS\n");

    for(int i=0;i<N;i++){
        SV* sv = svs1[i];
        pur(i) = (sv->prMes-sv->prInterp[sigInd])-pur0;
        printf("sv in N2:%d,%d,%f\n",sv->type,sv->svId,pur(i));

        Vector3d Iri = sv->position - referBase;
        Iri = Iri/Iri.norm();
        Vector3d Ir_0i = Ir0 - Iri;
        Gur(i,0) = Ir_0i(0);
        Gur(i,1) = Ir_0i(1);
        Gur(i,2) = Ir_0i(2);
        Gur(i,3) = double(sv->type==SYS_BDS);
    }
//    printf(" positon rtk unlock 0\n");
//    gnss->rtkManager.mtxData.unlock();
//    printf(" positon rtk unlock 1\n");
    cout<<"++++++++++++Gur\n"<<Gur<<endl;
    cout<<"++++++++++++pur\n"<<pur<<endl;


    MatrixXd GurT = Gur.transpose();
    Matrix4d H = (GurT*Gur).inverse();
//    cout<<"+++++H\n"<<H<<endl;

    Vector4d burt = H*(GurT*pur);
    cout<<"++++++bur\n"<<burt<<endl;

    xyz = referBase+burt.head(3);
    cout<<"++++++++++++basexyz\n"<<referBase<<endl;
    cout<<"++++++++++++xyz\n"<<xyz<<endl;

    XYZ2LLA(xyz,LLA);
    gnss->isPositioned = true;
    GNSS::PosRcd pos(rcvtow,xyz,vxyz);
    pos.PrintSol("+++++++xyz-end");

    fprintf(gnss->log,"xyz:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,xyz(0),xyz(1),xyz(2));
    fprintf(gnss->log,"LLA:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,LLA(1)*180/GPS_PI,LLA(0)*180/GPS_PI,LLA(2));

    double pc = (svCectre->position - xyz).norm();
    tu = svCectre->prMes+svCectre->tsDelta*Light_speed-svCectre->I-svCectre->T-pc;
    gnss->tu = tu;
    printf("tu= %f\n", tu);
    delete (this);
}

int PosSolver::PositionRtk() {
    int sigInd = 1;
    //插值
    int result = 0;

    vector<SV*> svs0[Nsys], svs1;
    PrepareSVsData(svs0);
    printf("rcvtow = %lf\n",rcvtow);
    //暂时单模计算
    if(numOfSys[SYS_BDS]*numOfSys[SYS_GPS]){
        printf("NOT one sys mode ------------count 0\n");
//        return -1;
    }
    SV* svCectre;
    int svCentreInd = -2;
    double elevMax = -5;
//    gnss->rtkManager.mtxData.lock();
    printf(" positon rtk lock 0\n");
    printf("---------calcu N0 = %d\n",svs0[7].size());

    for(int i = 0,k=-1; i<svs0[7].size();i++){
        SV* sv = svs0[7][i];
        double time = rcvtow - tu/Light_speed;
        if(sv->type==SYS_BDS)time-=14;
        double pr = sv->InterpRtkData(time,sigInd);
        printf("--sv:%d,%d-pr time= %f, interp = %f\n",
                sv->type,sv->svId,time,pr);

        if(0.0!=pr){
            svs1.push_back(sv);
            k++;
            printf("angle = %.2f\n",sv->elevationAngle);
            if(sv->elevationAngle>elevMax){
                elevMax = sv->elevationAngle;
                svCectre = sv;
                svCentreInd = k;
            }
        }
    }
    printf("---------calcu N1 = %d, svCentreInd = %d\n",svs1.size(),svCentreInd);
    if(svs1.size()<4){
        printf("N1 not enough");
        return -1;
    }

    svs1.erase(svs1.begin()+svCentreInd);

    int N = svs1.size();
    printf("---------calcu N1-1 = %d\n",N);
    if(N<3){
//        gnss->rtkManager.mtxData.unlock();
        return -1;
    }
    MatrixXd pur(N,1);
    MatrixXd Gur(N,3);

    double pur0 = svCectre->prMes - svCectre->prInterp[sigInd];
    printf("pur0-base:%f - %f\n",svCectre->prMes ,svCectre->prInterp[sigInd]);
    Vector3d referBase = gnss->rtkManager.ECEF_XYZ;
    Vector3d referLLA;
    XYZ2LLA(referBase,referLLA);
//    printf("++++++++++++BaseLLA === %lf,%lf,%lf\n",referLLA(1)*180/GPS_PI,referLLA(0)*180/GPS_PI,referLLA(2));

    Vector3d Ir0 = svCectre->position - referBase;
    Ir0 = Ir0/Ir0.norm();

    printf("start 2222222222222SSSSSSSSSSSSSSSSS\n");

    for(int i=0;i<N;i++){
        SV* sv = svs1[i];
        pur(i) = (sv->prMes-sv->prInterp[sigInd])-pur0;
        printf("sv in N2:%d,%d,%f\n",sv->type,sv->svId,pur(i));

        Vector3d svPositionEarthRotate;
        EarthRotate(sv->position,svPositionEarthRotate,(sv->position - xyz).norm()/Light_speed);

        Vector3d Iri = svPositionEarthRotate - referBase;
        Iri = Iri/Iri.norm();
        Vector3d Ir_0i = Ir0 - Iri;
        Gur(i,0) = Ir_0i(0);
        Gur(i,1) = Ir_0i(1);
        Gur(i,2) = Ir_0i(2);
    }
//    printf(" positon rtk unlock 0\n");
//    gnss->rtkManager.mtxData.unlock();
//    printf(" positon rtk unlock 1\n");
//    cout<<"++++++++++++Gur\n"<<Gur<<endl;
    cout<<"++++++++++++pur\n"<<pur<<endl;


    MatrixXd GurT = Gur.transpose();
    Matrix3d H = (GurT*Gur).inverse();
//    cout<<"+++++H\n"<<H<<endl;

    Vector3d bur = H*(GurT*pur);
    cout<<"++++++bur\n"<<bur<<endl;

    xyz = referBase+bur;
    cout<<"++++++++++++basexyz\n"<<referBase<<endl;
    cout<<"++++++++++++xyz\n"<<xyz<<endl;

    gnss->isPositioned = true;

    GNSS::PosRcd pos(rcvtow,xyz,vxyz);
    pos.PrintSol("+++++++xyz-end");
    gnss->AddPosRecord(pos);

    fprintf(gnss->log,"xyz:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,xyz(0),xyz(1),xyz(2));
    fprintf(gnss->log,"LLA:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,LLA(1)*180/GPS_PI,LLA(0)*180/GPS_PI,LLA(2));

    double pc = (svCectre->position - xyz).norm();
    tu = svCectre->prMes+svCectre->tsDelta*Light_speed-svCectre->I-svCectre->T-pc;
    gnss->tu = tu;
    printf("tu= %f\n", tu);
    char gga[128];
    MakeGGA(gga,LLA,rTime);
    rtk->SentGGA(gga,strlen(gga));
    delete (this);
}

int PosSolver::PositionSingle() {
    int result = 0;
    vector<SV*> svsForCalcu[Nsys];
    PrepareSVsData(svsForCalcu);

//    return 0;
    printf("rcvtow = %.4f\n",rcvtow);
    int N = svsForCalcu[7].size();
    if(N<4){
        printf("\n\n----------calcu:%d,Not enough Svs.\n",N);
        return -1;
    } else if(4==N && numOfSys[SYS_BDS]*numOfSys[SYS_GPS]){
        printf("4 SVs with SYS_bds and SYS_gps:%d, %d,   Unable to solve.\n",numOfSys[SYS_BDS],numOfSys[SYS_GPS]);
        return -1;
    } else if(0==numOfSys[SYS_BDS]*numOfSys[SYS_GPS]){
        result = SolvePosition(svsForCalcu[7]);

//        if(N>4){
//            for (int i = 0; i < N; ++i) {
//                vector<SV*> svsT = svsForCalcu;
//                svsT.erase(svsT.begin()+i);
//                SolvePosition(svsT);
//            }
//        }

    } else{
//        result = SolvePositionBeiDouGPS(svsForCalcu);
        result = SolvePosition(svsForCalcu[7]);
    }

//    if(result)
    delete(this);
}


int PosSolver::SolvePosition(vector<SV*>svsForCalcu) {
    int N = svsForCalcu.size();
    if(N<4)
        return -1;
    printf("\nSolve Position N = %d\n",N);
    Vector4d dtxyzt = Vector4d::Ones();
    MatrixXd pc(N,1);
    MatrixXd b(N,1);

    for(int i = 0;i< N;i++){
        SV *sv= svsForCalcu[i];
        pc(i) = sv->prMes + Light_speed * sv->tsDelta - sv->I - sv->T;
    }

    int numCalcu = 0;
    MatrixXd G(N,4);
    Matrix<double,4,4> H;
    while (dtxyzt.norm()>0.1){
//    while (numCalcu<20){
        numCalcu++;
        for(int i = 0;i< N;i++){
            SV *sv= svsForCalcu[i];
            double r = (sv->position-xyz).norm();
            //todo 自转修正
            Vector3d svPositionEarthRotate;
            EarthRotate(sv->position,svPositionEarthRotate,r/Light_speed);
            r = (svPositionEarthRotate-xyz).norm();

            b(i) = pc(i)-r-tu;
            G(i,0) = (xyz(0)-svPositionEarthRotate(0))/r;  //x
            G(i,1) = (xyz(1)-svPositionEarthRotate(1))/r;  //y
            G(i,2) = (xyz(2)-svPositionEarthRotate(2))/r;  //z
            G(i,3) = 1;
        }
        MatrixXd GT = G.transpose();
        H = (GT*G).inverse();
        dtxyzt = H*GT*b;

//        cout<<"b = "<<b<<endl;
        cout<<"dt = "<<dtxyzt<<endl;
        xyz += dtxyzt.head(3);
        tu+=dtxyzt(3);
//        cout<<"++++++++++++xyzt\n"<<xyz<<endl;
//        cout<<"xyzNorm\n"<<xyz.norm()<<endl;
//        XYZ2LLA(xyz,LLA);
//        cout<<"++++++++++++LLA\n"<<LLA<<endl;
        printf("num = %d, norm = %lf\n",numCalcu,dtxyzt.norm());
        if(numCalcu>30)return false;
    }
    double PDOP = sqrt(H(0,0)+H(1,1)+H(2,2));
    gnss->isPositioned = true;
    GNSS::PosRcd pos(rcvtow,xyz,vxyz);
    LLA = pos.lla;
    pos.PrintSol("+++++++xyz-From single one mode");
    gnss->AddPosRecord(pos);
    gnss->tu = tu;
//    fprintf(gnss->logDebug,"%f\n", tu);

    fprintf(gnss->log,"xyz:time = ,%.5f, pos = ,%.5f,%.5f,%.5f, PDOP=%.2f\n",rcvtow,xyz(0),xyz(1),xyz(2),PDOP);
    fprintf(gnss->log,"LLA:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,LLA(0)*180/GPS_PI,LLA(1)*180/GPS_PI,LLA(2));
    fprintf(gnss->log,"LLA:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,LLA(0),LLA(1),LLA(2));

//    cout<<"++++++++++++LLA\n"<<LLA*180/GPS_PI<<endl;

    printf("rcvtow = %lf\n",rcvtow);
    for(int i=0;i<N;i++)
    {
        double prres = b(i);
        SV *svTemp = svsForCalcu[i];
//        svTemp->CalcuelEvationAzimuth(xyz,LLA);
        printf("svsVis:time,%.4f,s%d%02d,pr,%.5f,prres,%.4f,%f\n",rcvtow,
               svTemp->type,svTemp->svId,svTemp->prMes,prres,svTemp->elevationAngle);
        fprintf(gnss->log,"svsVis:time,%.4f,s%d%02d,pr,%.5f,prres,%.4f,norm-x-y-z,%.5f,%.5f,%.5f,%.5f\n",rcvtow,
                svTemp->type,svTemp->svId,svTemp->prMes,prres,
                svTemp->position.norm(),svTemp->position(0),svTemp->position(1),svTemp->position(2));
    }
    
    printf("PDOP = %lf\n",PDOP);
    return 1;
}

int PosSolver::SolvePositionBeiDouGPS(vector<SV*>svsForCalcu){
    int N = svsForCalcu.size();
    printf("\nSolve Position SYS_GPS-SYS_BDS, N = %d\n",N);
    VectorXd dtxyzBG(5);
    Matrix<double, 5,5> H;
    dtxyzBG<<1,1,1,1,1;
    MatrixXd pc(N,1);
    MatrixXd b(N,1);

    for(int i = 0;i< N;i++){
        SV *sv= svsForCalcu[i];
        pc(i) = sv->prMes + Light_speed * sv->tsDelta - sv->I - sv->T;
    }

    int numCalcu = 0;
    while (dtxyzBG.norm()>0.1){
        numCalcu++;
        MatrixXd G(N,5);
        for(int i = 0;i< N;i++){
            SV *sv= svsForCalcu[i];
//            cout<<"rxyz=\n"<<xyz<<endl;
            double r = (sv->position-xyz).norm();
            //自转修正,这个我需要考虑一下是不是在这里处理
            double omegat = r/Light_speed*Omega_e;
            Matrix<double,3,3>earthRotate;
            earthRotate<<cos(omegat),sin(omegat),0,-sin(omegat),cos(omegat),0,0,0,1;
            Vector3d svPositionEarthRotate = earthRotate*sv->position;
//            cout<<"positon = \n"<<sv->position<<endl;
//            cout<<"earthRotate = \n"<<earthRotate<<endl;
//            cout<<"positionRotate = \n"<<svPositionEarthRotate<<endl;
            r = (svPositionEarthRotate-xyz).norm();
//            printf("r=%.10f\n",r);
//            r = (sv->position-rxyz).norm()+(sv->position(0)*rxyz(1)-sv->position(1)*rxyz(0))*Omega_e/Light_speed;
//            printf("r=%.10f\n",r);
//            printf("r=%.10f\n",(sv->position-rxyz).norm());

            G(i,0) = (xyz(0)-svPositionEarthRotate(0))/r;  //x
            G(i,1) = (xyz(1)-svPositionEarthRotate(1))/r;  //y
            G(i,2) = (xyz(2)-svPositionEarthRotate(2))/r;  //z
            G(i,3) = (SYS_BDS==sv->type)?1:0;
            G(i,4) = (SYS_GPS==sv->type)?1:0;
//            b(i) = pc(i)-r-tuBeiDou-tuGps;  //这里的r需要考虑一下
            b(i) = pc(i)-r-tuBeiDou*G(i,3)-tuGps*G(i,4);  //这里的r需要考虑一下

        }
        MatrixXd GT = G.transpose();
        H = (GT*G).inverse();
        dtxyzBG = H*GT*b;

//        cout<<"G=\n"<<G<<endl;
//        cout<<"GTG=\n"<<GT*G<<endl;
//        cout<<"GTGi=\n"<<(GT*G).inverse()<<endl;
//        cout<<"b = "<<b<<endl;
        cout<<"dt = "<<dtxyzBG<<endl;
        xyz += dtxyzBG.head(3);
        tuBeiDou+=dtxyzBG(3);
        tuGps+=dtxyzBG(4);
        printf("tuGps,tuBds %.15f,%.15f\n", tuGps / Light_speed, tuBeiDou / Light_speed);
        printf("num = %d, norm = %lf\n",numCalcu,dtxyzBG.norm());
        if(numCalcu>30)return false;
    }
    gnss->isPositioned = true;
    GNSS::PosRcd pos(rcvtow,xyz,vxyz);
    pos.PrintSol("+++++++xyz-end");
    gnss->AddPosRecord(pos);
    gnss->tuGps = tuGps;
    gnss->tuBeiDou = tuBeiDou;
    cout<<"++++++++++++rxyzt\n"<<xyz<<endl;
    printf("++++++++++++LLA === %lf,%lf,%lf\n",LLA(0)*180/GPS_PI,LLA(1)*180/GPS_PI,LLA(2));
    double PDOP = sqrt(H(0,0)+H(1,1)+H(2,2));
    printf("PDOP = %lf\n",PDOP);

    fprintf(gnss->log,"xyz:time = ,%.5f, pos = ,%.5f,%.5f,%.5f, PDOP=%.2f\n",rcvtow,xyz(0),xyz(1),xyz(2),PDOP);
    fprintf(gnss->log,"LLA:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,LLA(0)*180/GPS_PI,LLA(1)*180/GPS_PI,LLA(2));

    printf("rcvtow = %lf\n",rcvtow);
    for(int i=0;i<N;i++)
    {
        double prres = b(i);
        SV *svTemp = svsForCalcu[i];
//        svTemp->CalcuelEvationAzimuth(xyz,LLA);
        printf("svsVis:time,%.4f,s%d%02d,pr,%.5f,prres,%.4f,%f\n",rcvtow,
               svTemp->type,svTemp->svId,svTemp->prMes,prres,svTemp->elevationAngle);
        fprintf(gnss->log,"svsVis:time,%.4f,s%d%02d,pr,%.5f,prres,%.4f,norm-x-y-z,%.5f,%.5f,%.5f,%.5f\n",rcvtow,
                svTemp->type,svTemp->svId,svTemp->prMes,prres,
                svTemp->position.norm(),svTemp->position(0),svTemp->position(1),svTemp->position(2));
    }

    return 1;
}



int PosSolver::SolvePositionCalman() {

}

//todo move to serial data,and new class Measure;
int PosSolver::ReadVisibalSvsRaw(SVs *svs,vector<SV*> &svVisable, char *raw) {
    int maskGps[Ngps] = {0};
    int maskBds[Nbds] = {0};
    int tempTrack;
    char* playload = raw + 6;

    char* temp = playload;
    rcvtow = *(double*)playload;
    int week = *(uint16_t*)(playload+8);
    numMeas = *(u_int8_t*)(playload+11);
    rTime = GnssTime(week,rcvtow);
    printf("prepare rawdata , numMesa=%d\n",numMeas);
    if(0==numMeas)return -1;

    for(u_int8_t n = 0;n<numMeas;n++){
        SV* svTemp;
        int n32 = n*32;
        uint8_t gnssId = *(uint8_t *)(playload+36+n32);
        uint8_t svId = *(uint8_t *)(playload+37+n32);
        //watchable SVs
        switch (gnssId){
            case 0:
                svTemp = &(svs->svGpss[svId-1]);
                maskGps[svId-1] = 1;
                break;
            case 3:
                svTemp = &(svs->svBeiDous[svId-1]);
                maskBds[svId-1] = 1;
                break;
        }
        svVisable.push_back(svTemp);
        svTemp->prMes = *(double*)(playload+16+n32);
        svTemp->cpMes = *(double*)(playload+24+n32);
        svTemp->doMes = *(float *)(playload+32+n32);
//        printf("Measure______ %.5f,%.5f,%.5f\n", svTemp->prMes,svTemp->cpMes*19/100,svTemp->doMes);

        svTemp->measureRecord.push_back(SV::Measure(rcvtow,svTemp->prMes,svTemp->cpMes,svTemp->doMes));
    }
    for (int i = 0; i < Ngps; ++i) {
        tempTrack = svs->svGpss[i].trackingState+1;
        svs->svGpss[i].trackingState = maskGps[i]?tempTrack:0;
    }
    for (int i = 0; i < Nbds; ++i) {
        tempTrack = svs->svBeiDous[i].trackingState+1;
        svs->svBeiDous[i].trackingState = maskBds[i]?tempTrack:0;
    }
    return 0;
}

int PosSolver::SelectSvsFromVisible(vector<SV*> &all, vector<SV*> *select) {

    for(int i=0;i<numMeas;i++)
    {
        SV *svTemp = all[i];
        printf("getsv:%d,%02d:h:%d,pr=%lf,cp=%lf\t\t",svTemp->type,svTemp->svId,svTemp->SatH1,
                svTemp->prMes,svTemp->cpMes);

        //1,is used?
        if(!svTemp->open){
            printf("____closed\n");
            continue;
        }
        //2,judge ephemeric
        if(!svTemp->IsEphemOK(0,rTime))continue;
        if(!svTemp->IsEphemOK(gnss->ephemType,rTime))continue;
        //3,measure
        if(!svTemp->MeasureGood())continue;
        //4,elevtion angle
//        if(!svTemp->ElevGood())continue;
        //AODC?

//        if(svTemp->AODC>24){
//            printf("\n\n\nAAAAAODC = %d\n\n\n",svTemp->AODC);
//            useForCalcu = 0;
//        }
        //remain svs:
        select[7].push_back(svTemp);
        select[svTemp->type].push_back(svTemp);
        numOfSys[svTemp->type]++;
//        if(SV::SYS_BDS == svTemp->type){
//            numBDSUsed++;
//            select[1].push_back(svTemp);
//        }
//        if(SV::SYS_GPS == svTemp->type){
//            numGPSUsed++;
//            select[2].push_back(svTemp);
//        }
        printf("used this\n");
    }
}

int PosSolver::UpdateSvsPosition(vector<SV *> &svs, GnssTime rt, int ephType) {
    vector<SV *> result;
    printf("NNN== %d\n", svs.size());
    double ep[6];
    rt.time2epoch(ep);
    for(int i = 0;i< svs.size();i++){
        SV *sv= svs[i];

        //todo:rcvtow calculated from rt;
        double timeForECEF = rcvtow;
        if(SYS_BDS == sv->type)timeForECEF-=14;
        sv->CalcuTime(timeForECEF);
//        printf("ep = % 2.0f % 2.0f % 2.0f % 2.0f % 2.0f \n",ep[0],ep[1],ep[2],ep[3],ep[4]);
//        sv->PrintInfo(1);
//        printf("ts0 =  %.3f\n", sv->tsDelta*1e9);
//        sv->CalcuECEF(sv->tsReal);
//        printf("ts1 =  %.3f\n", sv->tsDelta*1e9);

        fprintf(gnss->log,"svsVisBrt:time,%.4f,s%d%02d,pr,%.5f,prres,%.4f,norm-x-y-z,%.5f,%.5f,%.5f,%.5f\n",rcvtow,
                sv->type,sv->svId,sv->prMes,0.0,
                sv->position.norm(),sv->position(0),sv->position(1),sv->position(2));

        Sp3Cell cell;
        GnssTime ts0 = rt;
        switch (ephType){
            case 1:
                ts0+=(-sv->prMes/Light_speed-sv->tsDelta);
                if(EphemSp3::Sp32ECEF(sv->ephemSp3->records, ts0, cell))break;
                printf("ts2 =  %.3f\n", cell.ts*1e9);

                sv->position = cell.pxyz;

                fprintf(gnss->log,"svsVisSp3:time,%.4f,s%d%02d,pr,%.5f,prres,%.4f,norm-x-y-z,%.5f,%.5f,%.5f,%.5f\n",rcvtow,
                        sv->type,sv->svId,sv->prMes,0.0,
                        sv->position.norm(),sv->position(0),sv->position(1),sv->position(2));
//                if(cell.ts!=0.0)sv->tsDelta = cell.ts;
//                sv->tsDelta-=sv->TGD1;
//                sv->PrintInfo(1);
                result.push_back(sv);
                break;
            case 0:
                if(!sv->IsEphemOK(0,rt))return -1;
                sv->CalcuECEF(sv->tsReal);
                result.push_back(sv);
//                sv->PrintInfo(1);
                break;
            default:
                continue;
        }

        if(gnss->isPositioned){
//            printf("elevation newed.\n");
            sv->CorrectIT(xyz,LLA,timeForECEF);
            sv->CalcuelEvationAzimuth(gnss->records.back().xyz,gnss->records.back().lla);
        } else{
//            printf("elevation default.\n");
            sv->CalcuelEvationAzimuth(gnss->xyzDefault,gnss->llaDefault);
        }

        XYZ2LLA(sv->position,sv->sLLA);
//        sv->PrintInfo(2);
//        sv->PrintInfo(1);
//        cout<<"norm"<<sv->position.norm()<<endl;
    }
    svs = result;
    return 0;
}

int PosSolver::MakeGGA(char *gga, Vector3d lla, GnssTime gpsTime) {
    GnssTime ggaTime = gpsTime;
    double h=0,ep[6],dms1[3],dms2[3],dop=1.0;
    int solq = 1;
    char *p=gga,*q,sum;

    ggaTime.gpst2utc();
    if (ggaTime.sec>=0.995) {ggaTime.time++; ggaTime.sec=0.0;}
    ggaTime.time2epoch(ep);

    deg2dms(fabs(lla(0))*180.0/GPS_PI,dms1);
    deg2dms(fabs(lla(1))*180.0/GPS_PI,dms2);
    p+=sprintf(p,"$GPGGA,%02.0f%02.0f%05.2f,%02.0f%010.7f,%s,%03.0f%010.7f,%s,%d,%02d,%.1f,%.3f,M,%.3f,M,%.1f,",
               ep[3],ep[4],ep[5],dms1[0],dms1[1]+dms1[2]/60.0,lla(0)>=0?"N":"S",
               dms2[0],dms2[1]+dms2[2]/60.0,lla(1)>=0?"E":"W",solq,
               6,dop,lla(2),h,1.0);
    for (q=(char *)gga+1,sum=0;*q;q++) sum^=*q; /* check-sum */
    p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
    return p-(char *)gga;
}

int PosSolver::PositionRtkKalman() {
    Vector3d xyzBackup = xyz;
    int sigInd = 1;
    int result = 0;
    Vector3d base = gnss->rtkManager.ECEF_XYZ;
    vector<SV*> svsForCalcu[Nsys];
    PrepareSVsData(svsForCalcu);

    SolvePosition(svsForCalcu[7]);
    printf("ttttttttttttttttttttttttttttttttu %f\n", tu);

    ProcessRtkData(svsForCalcu);
    Vector3d xyzSingle = xyz;

    int N = 0,xi=6,yi=0;
    for (int sys = 0; sys < Nsys-1; ++sys)N+=numOfSys[sys];
    printf("N= %d,nsys=%d\n", N,nsysUsed);
    if(N-nsysUsed<5){
        printf("svs not enough %d,nsys %d\n", N,nsysUsed);
        return -1;
    }
    VectorXd x(N+6),y(2*(N-nsysUsed)),hx(2*(N-nsysUsed));
    MatrixXd P(N+6,N+6),Q(N+6,N+6),Ppred(N+6,N+6);
    MatrixXd Hx(2*(N-nsysUsed),N+6),R(2*(N-nsysUsed),2*(N-nsysUsed));
    double cp_rb[2],pr_rb[2],r_rb[2],B_rb[2];
    vector<double*>pB,pPB,pRc,pRp;

    x.fill(0);y.fill(0);hx.fill(0);P.fill(0);Q.fill(0);Hx.fill(0);R.fill(0);
    x.head(6)<<xyzBackup,vxyz;
    P.block<6,6>(0,0)<<gnss->Pxv;
    Q.block<3,3>(0,0) = 1e5*Matrix3d::Identity();

    Vector3d bxlla,xyzrtk;
    //各系统
    for (int sys = 0,yHead,n,n_; sys < Nsys-1; ++sys) {
        n = numOfSys[sys];
        n_=n-1;
        if(0==n)continue;
        cout <<"sys="<< sys <<", n="<<n<<endl;

        double lambdai = GetFreq((SysType)sys, sigInd, true);
        printf("lambdai %f\n", lambdai);
        MatrixXd E(n,3),D(n-1,n),Rci(n,n),Rpi(n,n);
        D.fill(0);Rci.fill(0);Rpi.fill(0);
        //各卫星
        VectorXd yp(n-1);

        for (int i = 0,i_; i < n; ++i) {
            SV *sv = svsForCalcu[sys][i];
            Vector3d svPositionEarthRotate;
            EarthRotate(sv->position,svPositionEarthRotate,(sv->position - xyz).norm()/Light_speed);
            Vector3d e_sr = svPositionEarthRotate - xyz, e_sb = svPositionEarthRotate -base;
            E.block<1,3>(i,0)<<(e_sr/e_sr.norm()).transpose();

            printf("xi,yi %d,%d,\t", xi,yi);
            printf("rtksv:%d,%02d:h:%d,pr=%lf---%lf,.,cp=%lf,%.1f---%lf,%.1f [%.2f]\telev=%lf\n",sv->type,sv->svId,sv->SatH1,
                   sv->prMes,sv->prInterp[sigInd],sv->cpMes*lambdai,sv->cpMes-sv->prMes/lambdai,sv->cpInterp[sigInd],
                   (sv->cpInterp[sigInd]-sv->prInterp[sigInd])/lambdai,sv->cpMes-sv->cpInterp[sigInd]/lambdai,sv->elevationAngle);
            //单差部分
            double *tempBi=&(gnss->cycle[sys][sv->svId]);     pB.push_back(tempBi);
            double *tempPB=&(gnss->PB[sys][sv->svId]);        pPB.push_back(tempPB);
            double *tempRci=&(gnss->sigmaCy[sys][sv->svId]);  pRc.push_back(tempRci);
            double *tempRpi=&(gnss->sigmaPr[sys][sv->svId]);  pRp.push_back(tempRpi);

//            if(sv->KalmanFirst){
//                *tempBi = ((sv->cpMes*lambdai-sv->cpInterp[sigInd]) - (e_sr.norm()-e_sb.norm()) - tu)/lambdai;
//                sv->KalmanFirst = 0;
//            }

            double rdiff = sv->cpMes*lambdai-sv->prMes;
            double bdiff = sv->cpInterp[sigInd]-sv->prInterp[sigInd];
            double pr_ = sv->cpMes*lambdai-sv->cpInterp[sigInd];
            double pr__ = pr_ - tu;
            double r_ =  e_sr.norm()-e_sb.norm();
            double Biii = (rdiff-bdiff)/lambdai;
            printf("-----------cp- r-,tu,bi = %.2f,%.2f,%.2f,%.2f, %.2f\n",pr_,pr__,r_,tu,Biii);
            double rr = e_sr.norm(),rb = e_sb.norm();
            double ts = sv->tsDelta*Light_speed;
//            double tutemp = sv->prMes-rr+ts-sv->I-sv->T;
            double tub = sv->prInterp[sigInd]+ts-sv->I-sv->T-e_sb.norm();
            fprintf(gnss->logDebug,"%f,%f\n",tu,tub);
            double p_r_r = sv->prMes-rr-tu+ts;
            double p_r_b = sv->prInterp[sigInd] - rb - tub +ts;
//            fprintf(gnss->logDebug,"%.6f,\n",tu);
            if(sv->svId==0){
//                fprintf(gnss->logDebug,"%.6f,",sv->prInterp[sigInd]+sv->tsDelta*Light_speed-e_sb.norm());
//                fprintf(gnss->logDebug,"%.6f,",tu);
//                fprintf(gnss->logDebug,"%.6f,",sv->prMes+sv->tsDelta*Light_speed-e_sr.norm());
//
//                fprintf(gnss->logDebug,"%.6f,",Biii);
//                fprintf(gnss->logDebug,"%.6f,",sv->cpInterp[sigInd]);

                fprintf(gnss->logDebug,"%.6f,%.6f,",sv->I,sv->T);
                fprintf(gnss->logDebug,"%.6f,%.6f,",p_r_r,p_r_b);
                fprintf(gnss->logDebug,"%.6f,",tub);
                fprintf(gnss->logDebug,"%.6f,",tu);
                fprintf(gnss->logDebug,"%.6f,",rdiff);
                fprintf(gnss->logDebug,"%.6f,",bdiff);
                fprintf(gnss->logDebug,"%.6f,",(rdiff-bdiff)/lambdai);
                fprintf(gnss->logDebug,"\n");
            }
            //            *tempBi = Biii;
            x(xi) = *tempBi;
            P(xi,xi) = *tempPB;
            Rci(i,i) = 2*pow(*tempRci,2);
            Rpi(i,i) = 2*pow(*tempRpi,2);

            xi++;
            cp_rb[1] = sv->cpMes*lambdai-sv->cpInterp[sigInd];
            pr_rb[1] = sv->prMes-sv->prInterp[sigInd];
            B_rb[1] = *tempBi;
            r_rb[1] = e_sr.norm()-e_sb.norm();

            //双差部分
            if(0==i){
                cp_rb[0] = cp_rb[1];    pr_rb[0] = pr_rb[1];
                B_rb[0] = B_rb[1];      r_rb[0] = r_rb[1];

                yHead = 2*yi;
                continue;
            }
            i_ = i-1;
            D(i_,0)=1;D(i_,i)=-1;
            y(yHead+i_)     = cp_rb[0] - cp_rb[1];
            y(yHead+n_+i_)  = pr_rb[0] - pr_rb[1];
            yp(i_) = y(yHead+n_+i_);

            hx(yHead+i_)    = r_rb[0]-r_rb[1] + lambdai*(B_rb[0] - B_rb[1]);
            hx(yHead+n_+i_) = r_rb[0]-r_rb[1];

            printf("yi comp [%.2f,%.2f]\n", y(yHead+i_), y(yHead+n_+i_));
            printf("hxi comp [%.2f,%.2f]\n", hx(yHead+i_), hx(yHead+n_+i_));
            printf("B0,bi %f --- %f\n", B_rb[0],*tempBi);
            yi++;

            if(sv->svId==0){
                double rrb12 = r_rb[0]-r_rb[1];
                double prb12 = pr_rb[0]- pr_rb[1];
                fprintf(gnss->logDebug,"%f,%f,%f,%f\n",prb12,rrb12,prb12-rrb12,sv->tsDelta*Light_speed);
            }


        }
        Vector3d bx = (-D*E).colPivHouseholderQr().solve(yp);
        xyzrtk = bx+base;
//        x.head(3) = xyzrtk;
        GNSS::PosRcd pos00(rcvtow,xyzrtk,vxyz);
        pos00.PrintSol("+++++xyz 111111 ");

//        cout<<D<<E<<endl;
        MatrixXd Dt = D.transpose();
        Hx.block(yHead,0,n-1,3)<<-D*E;
        Hx.block(yHead+n-1,0,n-1,3)<<-D*E;
        Hx.block(yHead,xi-n,n-1,n)<<lambdai*D;
        R.block(yHead,yHead,n-1,n-1)<<D*Rci*Dt;
        R.block(yHead+n-1,yHead+n-1,n-1,n-1)<<D*Rpi*Dt;
    }
    P = 1*MatrixXd::Identity(N+6,N+6)+P;
//    P.block<3,3>(0,0) = P.block<3,3>(0,0)+0.1*Matrix3d::Identity();

    MatrixXd Ht = Hx.transpose();
    MatrixXd Kk = P*Ht*((Hx*P*Ht+R).inverse());
    Ppred = (MatrixXd::Identity(N+6,N+6)-Kk*Hx)*P;
    VectorXd xyzPred = x + Kk*(y-hx);

    cout <<"\nx\n"<<x.transpose() <<"\ny\n"<<y.transpose()<<"\nhx\n"<<hx.transpose()<<"\ny-hx\n"<<(y-hx).transpose()<<endl;
//    cout<<"\nP\n"<<P <<"\nHx\n"<<Hx<<"\npPred\n"<<Ppred<<endl;
//    cout <<"\nKk\n"<<Kk<<endl;
//    cout <<"\nR\n"<<R<<endl<<;
    cout<<"\nxyzAdd\n"<<(Kk*(y-hx)).transpose() <<endl;
    cout<<"\nxyzPred\n"<<xyzPred.transpose() <<endl;
    //更新
    xyz = xyzPred.head(3);
//    xyz = xyzrtk;
    vxyz = xyzPred.segment(3,3);
    for (int i = 0; i < pB.size(); ++i) {
        *pB[i] = xyzPred(6+i);
        *pPB[i] = Ppred(6+i,6+i);
    }
    gnss->Pxv = Ppred.block<6,6>(0,0);
//    gnss->Pxv(0,0) = Ppred(0,0);
//    gnss->Pxv(1,1) = Ppred(1,1);
//    gnss->Pxv(2,2) = Ppred(2,2);

    //显示
    gnss->isPositioned = true;
    GNSS::PosRcd pos(rcvtow,xyz,vxyz);
    pos.PrintSol("+++++++xyz-222222");
    gnss->AddPosRecord(pos);
    cout<<"++++xyzDiff single\n"<<-gnss->xyz00<<endl;
    cout<<"++++xyzDiff rtk\n"<<xyzrtk-gnss->xyz00<<endl;
    cout<<"++++xyzDiff rtk kalman\n"<<xyzPred.head(3)-gnss->xyz00<<endl;

    fprintf(gnss->log,"xyz:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,xyz(0),xyz(1),xyz(2));
    fprintf(gnss->log,"LLA:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,LLA(1)*180/GPS_PI,LLA(0)*180/GPS_PI,LLA(2));

}

int PosSolver::ProcessRtkData(vector<SV *> *select) {
    int sigInd = 1;
    nsysUsed = 0;
    SV* svCectre;
    select[7].clear();
    for(int sys = 0;sys<Nsys-1;sys++){
        if(select[sys].empty())continue;
        double elevMax = -5;
        int svCentreInd = -2;
        int nAll = select[sys].size();
        vector<SV*>svsTemp;
        int nRtk = 0;
        for (int i = 0; i < nAll; ++i) {
            SV *sv = select[sys][i];
            double time = rcvtow - tu/Light_speed;
            if(sv->type==SYS_BDS)time-=14;

            double pr = sv->InterpRtkData(time,sigInd);
            if(0!=pr){
                nRtk++;
                svsTemp.push_back(sv);
                if(sv->elevationAngle>elevMax){
                    elevMax = sv->elevationAngle;
                    svCectre = sv;
                    svCentreInd = nRtk-1;
                }
            }
        }
        if(!svsTemp.empty()){
            svsTemp.erase(svsTemp.begin()+svCentreInd);
            svsTemp.insert(svsTemp.begin(),svCectre);

//            Vector3d xyzrot;
//            double pc = (svCectre->position - xyz).norm();
//            EarthRotate(svCectre->position,xyzrot,pc/Light_speed);
//            pc = (xyzrot-xyz).norm();
//            tu = svCectre->prMes+svCectre->tsDelta*Light_speed-svCectre->I-svCectre->T-pc;
////            fprintf(gnss->logDebug,"%f\n",tu);
//            gnss->tu = tu;
        }
        select[sys].swap(svsTemp);
        numOfSys[sys]=0;
        if(nRtk>1){
            numOfSys[sys] = nRtk;
            nsysUsed++;
        }
    }

}