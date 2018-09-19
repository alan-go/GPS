#include "PosSolver.h"
#include "GNSS.h"
#include "EphemSp3.h"

PosSolver::PosSolver(){}

PosSolver::PosSolver(SVs *svs, NtripRTK *rtk, GNSS *gnss): svs(svs),rtk(rtk),gnss(gnss),numBDSUsed(0),numGPSUsed(0){
    if(gnss->isPositioned){
        xyz = gnss->records[0].xyz;
        LLA = gnss->records[0].lla;
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

    numGPSUsed = 0;numBDSUsed = 0;
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
    UpdateSvsPosition(svsOut[0], rTime, gnss->ephemType);

}

int PosSolver::PositionRtk2() {
    int sigInd = 1;
    //插值
    int result = 0;

    vector<SV*> svs0[3], svs1;
    PrepareSVsData(svs0);
    printf("rcvtow = %lf\n",rcvtow);
    //暂时单模计算
    if(numGPSUsed*numBDSUsed){
        printf("NOT one sys mode ------------count 0\n");
//        return -1;
    }
    SV* svCectre;
    int svCentreInd = -2;
    double elevMax = -5;
//    gnss->rtkManager.mtxData.lock();
    printf(" positon rtk lock 0\n");
    printf("---------calcu N0 = %d\n",svs0[0].size());

    for(int i = 0,k=-1; i<svs0[0].size();i++){
        SV* sv = svs0[0][i];
        double time = rcvtow - tu/Light_speed;
        if(sv->type==SV::SYS_BDS)time-=14;
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
        Gur(i,3) = double(sv->type==SV::SYS_BDS);
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
    gnss->AddPosRecord(GNSS::PosRcd(rcvtow,xyz,LLA));
    printf("++++++++++++LLA === %lf,%lf,%lf\n",LLA(1)*180/GPS_PI,LLA(0)*180/GPS_PI,LLA(2));

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

    vector<SV*> svs0[3], svs1;
    PrepareSVsData(svs0);
    printf("rcvtow = %lf\n",rcvtow);
    //暂时单模计算
    if(numGPSUsed*numBDSUsed){
        printf("NOT one sys mode ------------count 0\n");
//        return -1;
    }
    SV* svCectre;
    int svCentreInd = -2;
    double elevMax = -5;
//    gnss->rtkManager.mtxData.lock();
    printf(" positon rtk lock 0\n");
    printf("---------calcu N0 = %d\n",svs0[0].size());

    for(int i = 0,k=-1; i<svs0[0].size();i++){
        SV* sv = svs0[0][i];
        double time = rcvtow - tu/Light_speed;
        if(sv->type==SV::SYS_BDS)time-=14;
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

        Vector3d Iri = sv->position - referBase;
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

    XYZ2LLA(xyz,LLA);
    gnss->isPositioned = true;
    gnss->AddPosRecord(GNSS::PosRcd(rcvtow,xyz,LLA));
    printf("++++++++++++LLA === %lf,%lf,%lf\n",LLA(1)*180/GPS_PI,LLA(0)*180/GPS_PI,LLA(2));

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
    vector<SV*> svsForCalcu[3];
    PrepareSVsData(svsForCalcu);

//    return 0;
    printf("rcvtow = %.4f\n",rcvtow);
    int N = svsForCalcu[0].size();
    if(N<4){
        printf("\n\n----------calcu:%d,Not enough Svs.\n",N);
        return -1;
    } else if(4==N && numGPSUsed*numBDSUsed){
        printf("4 SVs with SYS_GPS and SYS_BDS:%d, %d,   Unable to solve.\n",numGPSUsed,numBDSUsed);
        return -1;
    } else if(0==numBDSUsed*numGPSUsed){
        result = SolvePosition(svsForCalcu[0]);

//        if(N>4){
//            for (int i = 0; i < N; ++i) {
//                vector<SV*> svsT = svsForCalcu;
//                svsT.erase(svsT.begin()+i);
//                SolvePosition(svsT);
//            }
//        }

    } else{
//        result = SolvePositionBeiDouGPS(svsForCalcu);
        result = SolvePosition(svsForCalcu[0]);
    }

//    if(result)
    delete(this);
}


int PosSolver::SolvePosition(vector<SV*>svsForCalcu) {
    int N = svsForCalcu.size();
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
    XYZ2LLA(xyz,LLA);
    gnss->isPositioned = true;
    gnss->AddPosRecord(GNSS::PosRcd(rcvtow,xyz,LLA));
    gnss->tu = tu;
    printf("tu= %f\n", tu);
    cout<<"++++++++++++xyz\n"<<xyz<<endl;
    printf("++++++++++++LLA === %lf,%lf,%lf\n",LLA(0)*180/GPS_PI,LLA(1)*180/GPS_PI,LLA(2));

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
            G(i,3) = (SV::SYS_BDS==sv->type)?1:0;
            G(i,4) = (SV::SYS_GPS==sv->type)?1:0;
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
    XYZ2LLA(xyz,LLA);
    gnss->isPositioned = true;
    gnss->AddPosRecord(GNSS::PosRcd(rcvtow,xyz,LLA));
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
    int maskGps[NGPS] = {0};
    int maskBds[NBeiDou] = {0};
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
    for (int i = 0; i < NGPS; ++i) {
        tempTrack = svs->svGpss[i].trackingState+1;
        svs->svGpss[i].trackingState = maskGps[i]?tempTrack:0;
    }
    for (int i = 0; i < NBeiDou; ++i) {
        tempTrack = svs->svBeiDous[i].trackingState+1;
        svs->svBeiDous[i].trackingState = maskBds[i]?tempTrack:0;
    }
    return 0;
}

int PosSolver::SelectSvsFromVisible(vector<SV*> &all, vector<SV*> *select) {

    for(int i=0;i<numMeas;i++)
    {
        SV *svTemp = all[i];
        printf("svs visiable:%d,%02d:%d,%d,%d.health:%d,pr=%lf\n",
               svTemp->type,svTemp->svId,svTemp->bstEphemOK[0],svTemp->bstEphemOK[1],
               svTemp->bstEphemOK[2],svTemp->SatH1,svTemp->prMes);

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
        select[0].push_back(svTemp);
//        numBDSUsed+=(SV::SYS_BDS == svTemp->type);
//        numGPSUsed+=(SV::SYS_GPS == svTemp->type);
        if(SV::SYS_BDS == svTemp->type){
            numBDSUsed++;
            select[1].push_back(svTemp);
        }
        if(SV::SYS_GPS == svTemp->type){
            numGPSUsed++;
            select[2].push_back(svTemp);
        }
        printf("svsfor calcu\n");
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
        if(SV::SYS_BDS == sv->type)timeForECEF-=14;
        sv->CalcuTime(timeForECEF);
        printf("ep = % 2.0f % 2.0f % 2.0f % 2.0f % 2.0f \n",ep[0],ep[1],ep[2],ep[3],ep[4]);
        sv->PrintInfo(1);
        printf("ts0 =  %.3f\n", sv->tsDelta*1e9);
        sv->CalcuECEF(sv->tsReal);
        printf("ts1 =  %.3f\n", sv->tsDelta*1e9);

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
            printf("elevation newed.\n");
            sv->CorrectIT(xyz,LLA,timeForECEF);
            sv->CalcuelEvationAzimuth(gnss->records[0].xyz,gnss->records[0].lla);
        } else{
            printf("elevation default.\n");
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
    int result = 0;
    vector<SV*> svsForCalcu[3];
    PrepareSVsData(svsForCalcu);
    SelectRtkData(svsForCalcu);
}

int PosSolver::SelectRtkData(vector<SV *> *select) {
    int sigInd = 1;
    SV* svCectre;

    for(int sys = 1;sys<3;sys++){
        double elevMax = -5;
        int svCentreInd = -2;
        int nsys = select[sys].size();
        vector<SV*>svsTemp;
        int k = -1;
        for (int i = 0; i < nsys; ++i) {
            SV *sv = select[sys][i];
            double time = rcvtow - tu/Light_speed;
            if(sv->type==SV::SYS_BDS)time-=14;
            double pr = sv->InterpRtkData(time,sigInd);
            if(pr){
                k++;
                svsTemp.push_back(sv);
                if(sv->elevationAngle>elevMax){
                    elevMax = sv->elevationAngle;
                    svCectre = sv;
                    svCentreInd = k;
                }
            }
        }
        if(!svsTemp.empty()){
            svsTemp.erase(svsTemp.begin()+svCentreInd);
            svsTemp.insert(svsTemp.begin(),svCectre);
        }
        select[sys].swap(svsTemp);
    }
    select[0].clear();
}