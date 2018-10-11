#include "PosSolver.h"
#include "GNSS.h"
#include "EphemSp3.h"

auto Diff_Vec = [](VectorXd x,int n)->VectorXd { return x(0)*VectorXd::Ones(n-1)-x.tail(n-1);};

auto Diff_Pr_rb = [](SV* sv,int sigId)->double { return sv->measureDat.front()->prMes-sv->prInterp[sigId];};
auto Diff_Cp_rb = [](SV* sv,int sigId)->double { return sv->measureDat.front()->cpMes-sv->cpInterp[sigId];};
auto Diff_Cyc_0 = [](SV* sv,int sigId)->double { return sv->measureDat.front()->cycle;};
auto Diff_cyP_0 = [](SV* sv,int sigId)->double { return sv->measureDat.front()->cycleP;};
auto Diff_Rcp_0 = [](SV* sv,int sigId)->double { return sv->measureDat.front()->sigmaCp;};
auto Diff_Rpr_0 = [](SV* sv,int sigId)->double { return sv->measureDat.front()->sigmaPr;};


PosSolver::PosSolver(){}

PosSolver::PosSolver(SvAll svs, NtripRTK *rtk, GNSS *gnss):svsBox(svs),rtk(rtk),gnss(gnss){
    Solution recent = gnss->records.front();
    if(gnss->isPositioned){
        xyz = recent.xyz;
//        LLA = gnss->records.back().lla;
    } else{
        xyz = gnss->xyzDefault;
//        LLA = gnss->llaDefault;
    }
    for (int i = 0; i < Nsys; ++i) {
        tu[i] = recent.tu[i];
    }
}

PosSolver::~PosSolver(){
    for(SvSys* sys:svsBox.sysUsed)delete(sys);
}

int PosSolver::PrepareSVsData(vector<SV*> &svsIn) {
    SelectSvsFromVisible(svsIn);
    if(svsBox.svUsedAll.empty()){ return -1;}
    timeSolver = svsBox.svUsedAll.front()->measureDat.front()->time;
    UpdateSvsPosition(svsBox.svUsedAll, timeSolver, gnss->ephemType);

}


int PosSolver::PositionSingle(vector<SV*> _svsIn) {
    PrepareSVsData(_svsIn);
    printf("single sinlel========== count=%d,nSat=%d,nSys=%d\n", gnss->count,nSat,nSys);

    soltion = Solution(timeSolver,xyz,vxyz);
    soltion.Show("0000000before");
    if(nSat-nSys<4){
        printf("Not enough svs nSat,nSys=%d,%d\n",nSat,nSys);
        return -1;
    }
    double threshold = 0.1;
    VectorXd dt_xyz_tu(3+nSys),pc(nSat),b(nSat);
    MatrixXd H(3+nSys,3+nSys);
    dt_xyz_tu.fill(1);
    int i =0,count =0;

//    for(SV* sv:svsBox.svUsedAll)
//        pc(i++) = sv->measureDat.front()->prMes + Light_speed * sv->tsDelta - sv->I - sv->T;

    while (dt_xyz_tu.norm()>threshold){
        int yhead=0,xhead=3;
        MatrixXd G(nSat,3+nSys),Gt(3+nSys,nSat);
        G.fill(0);
        for(SvSys* sys:svsBox.sysUsed){
            i=0;
            int ntemp = sys->table.size();
            if(0==ntemp)continue;
            MatrixXd E(ntemp,3);
            for(SV* sv:sys->table){
                double tRotate = ((sv->position-xyz).norm()+sv->I+sv->T)/Light_speed;
                Vector3d posRotate,e;
                EarthRotate(sv->position,posRotate,tRotate);
                e = posRotate-xyz;
                double r = e.norm();
                E.block<1,3>(i,0)<<(e/r).transpose();
//                printf("yhead,i %d,%d\n", yhead, i);
//                sv->PrintInfo(1);
                double pc = sv->measureDat.front()->prMes + Light_speed * sv->tsDelta - sv->I - sv->T;
                b(yhead+i)=pc-r-tu[sys->type];
//                printf("b---- %.2f,%.2f,%.2f\n",pc,r,tu[sys->type] );
                i++;
            }
            G.block(yhead,0,ntemp,3)<<-E;
            G.block(yhead,xhead,ntemp,1)<<VectorXd::Ones(ntemp);
            yhead+=ntemp;xhead++;
        }
//        dt_xyz_tu =  G.colPivHouseholderQr().solve(b);
        Gt = G.transpose();
        H = (Gt*G).inverse();
        dt_xyz_tu = H*Gt*b;
        xyz+=dt_xyz_tu.head(3);
//        cout<<"dt = "<<endl<<dt_xyz_tu<<endl;
//        cout<<"b = "<<endl<<b<<endl;
//        cout<<"G= "<<endl<<G<<endl;
//        cout<<"pc = "<<endl<<pc<<endl;

        i=3;
        //todo:sysused:
        for(SvSys* sys:svsBox.sysUsed)
            if(!sys->table.empty())
                tu[sys->type]+=dt_xyz_tu(i++);
        if(++count>10)return -1;
    }
    PDOP = sqrt(H(0,0)+H(1,1)+H(2,2));
    soltion = Solution(timeSolver,xyz,vxyz);
    return 0;
}



int PosSolver::SelectSvsFromVisible(vector<SV*> &all) {
    for(SV* sv:all)
    {
        Measure *msTemp = sv->measureDat.front();
        printf("getsv:%d,%02d:h:%d,pr=%lf,cp=%lf\t\t",sv->type,sv->svId,sv->SatH1,
                msTemp->prMes,msTemp->cpMes);

        //1,is used?
        if(!sv->open){
            printf("____closed\n");
            continue;
        }
        //2,judge ephemeric
        if(!sv->IsEphemOK(0,timeSolver))continue;
        if(!sv->IsEphemOK(gnss->ephemType,timeSolver))continue;
        //3,measure
        if(!sv->MeasureGood())continue;
        //4,elevtion angle
//        if(!sv->ElevGood())continue;
        //AODC?

//        if(sv->AODC>24){
//            printf("\n\n\nAAAAAODC = %d\n\n\n",sv->AODC);
//            useForCalcu = 0;
//        }
        //remain svs:


        svsBox.AddToUsed(sv);
       printf("+++this one useable\n");
    }
    nSat = svsBox.svUsedAll.size();
    nSys = svsBox.sysUsed.size();
}

int PosSolver::UpdateSvsPosition(vector<SV *> &svs, GnssTime rt, int ephType) {
    vector<SV *> result;
    printf("N0= %d\n", svs.size());
    printf("rcvtow %.3f\n", rt.tow);
    double ep[6];
    rt.time2epoch(ep);
    for(SV* sv:svs){
        Measure *measure = sv->measureDat.front();
        //todo:rcvtow calculated from rt;
        double timeForECEF = rt.tow;
        if(SYS_BDS == sv->type)timeForECEF-=14;
        sv->CalcuTime(timeForECEF);
//        printf("ep = % 2.0f % 2.0f % 2.0f % 2.0f % 2.0f \n",ep[0],ep[1],ep[2],ep[3],ep[4]);
//        sv->PrintInfo(1);
//        printf("ts0 =  %.3f\n", sv->tsDelta*1e9);
//        sv->CalcuECEF(sv->tsReal);
//        printf("ts1 =  %.3f\n", sv->tsDelta*1e9);

//        fprintf(gnss->log,"svsVisBrt:time,%.4f,s%d%02d,pr,%.5f,prres,%.4f,norm-x-y-z,%.5f,%.5f,%.5f,%.5f\n",rcvtow,
//                sv->type,sv->svId,measure->prMes,0.0,
//                sv->position.norm(),sv->position(0),sv->position(1),sv->position(2));

        Sp3Cell cell;
        GnssTime ts0 = rt;
        switch (ephType){
            case 1:
                ts0+=(-measure->prMes/Light_speed-sv->tsDelta);
                if(EphemSp3::Sp32ECEF(sv->ephemSp3->records, ts0, cell))break;
                printf("ts2 =  %.3f\n", cell.ts*1e9);

                sv->position = cell.pxyz;

//                fprintf(gnss->log,"svsVisSp3:time,%.4f,s%d%02d,pr,%.5f,prres,%.4f,norm-x-y-z,%.5f,%.5f,%.5f,%.5f\n",rcvtow,
//                        sv->type,sv->svId,measure->prMes,0.0,
//                        sv->position.norm(),sv->position(0),sv->position(1),sv->position(2));
//                if(cell.ts!=0.0)sv->tsDelta = cell.ts;
//                sv->tsDelta-=sv->TGD1;
//                sv->PrintInfo(1);
                result.push_back(sv);
                break;
            case 0:
                if(!sv->IsEphemOK(0,rt))return -1;
                sv->CalcuECEF(sv->tsReal);
//                sv->PrintInfo(0);
//                cout<<"svpositon"<<endl<<sv->position<<endl;
                result.push_back(sv);
//                sv->PrintInfo(1);
                break;
            default:
                continue;
        }

        if(gnss->isPositioned){
//            printf("elevation newed.\n");
            sv->CorrectIT(xyz,timeForECEF);
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



int PosSolver::ProcessRtkData() {
    int sigInd = 1;
    nSat=nSys=0;
    SvAll tempAll = gnss->svsManager;
    svsBox.svUsedAll.clear();

    for(SvSys* sys:svsBox.sysUsed){
        double elevMax = -5;
        int svCentreInd = -2;
        SV* svCectre;
        vector<SV*> temp,temp0;
        temp.swap(sys->table);
        int n = temp.size(),nRtk = 0;
        for(SV* sv:temp){
            double time = timeSolver.tow - tu[sv->type]/Light_speed;
            if(sv->type==SYS_BDS)time-=14;

            double pr = sv->InterpRtkData(time,sigInd);
            if(0!=pr){
                nRtk++;
                temp0.push_back(sv);
                if(sv->elevationAngle>elevMax){
                    elevMax = sv->elevationAngle;
                    svCectre = sv;
                    svCentreInd = nRtk-1;
                }
            }
        }
        if(temp0.size()>1){
            temp0.erase(temp0.begin()+svCentreInd);
            svsBox.AddToUsed(svCectre);
            for(SV* svv:temp0)svsBox.AddToUsed(svv);
        }
    }
    nSat = svsBox.svUsedAll.size();
    nSys = svsBox.sysUsed.size();
}

//int PosSolver::PositionRtk2() {
//    int sigInd = 1;
//    //插值
//    int result = 0;
//
//    vector<SV*> svs0[Nsys], svs1;
//    PrepareSVsData(svs0);
//    printf("rcvtow = %lf\n",rcvtow);
//    //暂时单模计算
//    if(sysEachNum[SYS_BDS]*sysEachNum[SYS_GPS]){
//        printf("NOT one sys mode ------------count 0\n");
////        return -1;
//    }
//    SV* svCectre;
//    int svCentreInd = -2;
//    double elevMax = -5;
////    gnss->rtkManager.mtxData.lock();
//    printf(" positon rtk lock 0\n");
//    printf("---------calcu N0 = %d\n",svs0[7].size());
//
//    for(int i = 0,k=-1; i<svs0[7].size();i++){
//        SV* sv = svs0[7][i];
//        double time = rcvtow - tui/Light_speed;
//        if(sv->type==SYS_BDS)time-=14;
//        double pr = sv->InterpRtkData(time,sigInd);
//        printf("--sv:%d,%d-pr time= %f, interp = %f\n",
//               sv->type,sv->svId,time,pr);
//
//        if(0.0!=pr){
//            svs1.push_back(sv);
//            k++;
//            printf("angle = %.2f\n",sv->elevationAngle);
//            if(sv->elevationAngle>elevMax){
//                elevMax = sv->elevationAngle;
//                svCectre = sv;
//                svCentreInd = k;
//            }
//        }
//    }
//    printf("---------calcu N1 = %d, svCentreInd = %d\n",svs1.size(),svCentreInd);
//    if(svs1.size()<4){
//        printf("N1 not enough");
//        return -1;
//    }
//
//    svs1.erase(svs1.begin()+svCentreInd);
//
//    int N = svs1.size();
//    printf("---------calcu N2 = %d\n",N);
//    if(N<4){
////        gnss->rtkManager.mtxData.unlock();
//        return -1;
//    }
//    MatrixXd pur(N,1);
//    MatrixXd Gur(N,4);
//
//    double pur0 = svCectre->prMes - svCectre->prInterp[sigInd];
//    printf("pur0-base:%f - %f\n",svCectre->prMes ,svCectre->prInterp[sigInd]);
//    Vector3d referBase = gnss->rtkManager.ECEF_XYZ;
//    Vector3d referLLA;
//    XYZ2LLA(referBase,referLLA);
////    printf("++++++++++++BaseLLA === %lf,%lf,%lf\n",referLLA(1)*180/GPS_PI,referLLA(0)*180/GPS_PI,referLLA(2));
//
//    Vector3d Ir0 = svCectre->position - referBase;
//    Ir0 = Ir0/Ir0.norm();
//
//    printf("start 2222222222222SSSSSSSSSSSSSSSSS\n");
//
//    for(int i=0;i<N;i++){
//        SV* sv = svs1[i];
//        pur(i) = (sv->prMes-sv->prInterp[sigInd])-pur0;
//        printf("sv in N2:%d,%d,%f\n",sv->type,sv->svId,pur(i));
//
//        Vector3d Iri = sv->position - referBase;
//        Iri = Iri/Iri.norm();
//        Vector3d Ir_0i = Ir0 - Iri;
//        Gur(i,0) = Ir_0i(0);
//        Gur(i,1) = Ir_0i(1);
//        Gur(i,2) = Ir_0i(2);
//        Gur(i,3) = double(sv->type==SYS_BDS);
//    }
////    printf(" positon rtk unlock 0\n");
////    gnss->rtkManager.mtxData.unlock();
////    printf(" positon rtk unlock 1\n");
//    cout<<"++++++++++++Gur\n"<<Gur<<endl;
//    cout<<"++++++++++++pur\n"<<pur<<endl;
//
//
//    MatrixXd GurT = Gur.transpose();
//    Matrix4d H = (GurT*Gur).inverse();
////    cout<<"+++++H\n"<<H<<endl;
//
//    Vector4d burt = H*(GurT*pur);
//    cout<<"++++++bur\n"<<burt<<endl;
//
//    xyz = referBase+burt.head(3);
//    cout<<"++++++++++++basexyz\n"<<referBase<<endl;
//    cout<<"++++++++++++xyz\n"<<xyz<<endl;
//
//    XYZ2LLA(xyz,LLA);
//    gnss->isPositioned = true;
//    Solution pos(rcvtow,xyz,vxyz);
//    pos.Show("+++++++xyz-end");
//
//    fprintf(gnss->log,"xyz:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,xyz(0),xyz(1),xyz(2));
//    fprintf(gnss->log,"LLA:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,LLA(1)*180/GPS_PI,LLA(0)*180/GPS_PI,LLA(2));
//
//    double pc = (svCectre->position - xyz).norm();
//    tui = svCectre->prMes+svCectre->tsDelta*Light_speed-svCectre->I-svCectre->T-pc;
//    gnss->tu = tui;
//    printf("tui= %f\n", tui);
//    delete (this);
//}
//
//int PosSolver::PositionRtk() {
//    int sigInd = 1;
//    //插值
//    int result = 0;
//
//    vector<SV*> svs0[Nsys], svs1;
//    PrepareSVsData(svs0);
//    printf("rcvtow = %lf\n",rcvtow);
//    //暂时单模计算
//    if(sysEachNum[SYS_BDS]*sysEachNum[SYS_GPS]){
//        printf("NOT one sys mode ------------count 0\n");
////        return -1;
//    }
//    SV* svCectre;
//    int svCentreInd = -2;
//    double elevMax = -5;
////    gnss->rtkManager.mtxData.lock();
//    printf(" positon rtk lock 0\n");
//    printf("---------calcu N0 = %d\n",svs0[7].size());
//
//    for(int i = 0,k=-1; i<svs0[7].size();i++){
//        SV* sv = svs0[7][i];
//        double time = rcvtow - tui/Light_speed;
//        if(sv->type==SYS_BDS)time-=14;
//        double pr = sv->InterpRtkData(time,sigInd);
//        printf("--sv:%d,%d-pr time= %f, interp = %f\n",
//                sv->type,sv->svId,time,pr);
//
//        if(0.0!=pr){
//            svs1.push_back(sv);
//            k++;
//            printf("angle = %.2f\n",sv->elevationAngle);
//            if(sv->elevationAngle>elevMax){
//                elevMax = sv->elevationAngle;
//                svCectre = sv;
//                svCentreInd = k;
//            }
//        }
//    }
//    printf("---------calcu N1 = %d, svCentreInd = %d\n",svs1.size(),svCentreInd);
//    if(svs1.size()<4){
//        printf("N1 not enough");
//        return -1;
//    }
//
//    svs1.erase(svs1.begin()+svCentreInd);
//
//    int N = svs1.size();
//    printf("---------calcu N1-1 = %d\n",N);
//    if(N<3){
////        gnss->rtkManager.mtxData.unlock();
//        return -1;
//    }
//    MatrixXd pur(N,1);
//    MatrixXd Gur(N,3);
//
//    double pur0 = svCectre->prMes - svCectre->prInterp[sigInd];
//    printf("pur0-base:%f - %f\n",svCectre->prMes ,svCectre->prInterp[sigInd]);
//    Vector3d referBase = gnss->rtkManager.ECEF_XYZ;
//    Vector3d referLLA;
//    XYZ2LLA(referBase,referLLA);
////    printf("++++++++++++BaseLLA === %lf,%lf,%lf\n",referLLA(1)*180/GPS_PI,referLLA(0)*180/GPS_PI,referLLA(2));
//
//    Vector3d Ir0 = svCectre->position - referBase;
//    Ir0 = Ir0/Ir0.norm();
//
//    printf("start 2222222222222SSSSSSSSSSSSSSSSS\n");
//
//    for(int i=0;i<N;i++){
//        SV* sv = svs1[i];
//        pur(i) = (sv->prMes-sv->prInterp[sigInd])-pur0;
//        printf("sv in N2:%d,%d,%f\n",sv->type,sv->svId,pur(i));
//
//        Vector3d svPositionEarthRotate;
//        EarthRotate(sv->position,svPositionEarthRotate,(sv->position - xyz).norm()/Light_speed);
//
//        Vector3d Iri = svPositionEarthRotate - referBase;
//        Iri = Iri/Iri.norm();
//        Vector3d Ir_0i = Ir0 - Iri;
//        Gur(i,0) = Ir_0i(0);
//        Gur(i,1) = Ir_0i(1);
//        Gur(i,2) = Ir_0i(2);
//    }
////    printf(" positon rtk unlock 0\n");
////    gnss->rtkManager.mtxData.unlock();
////    printf(" positon rtk unlock 1\n");
////    cout<<"++++++++++++Gur\n"<<Gur<<endl;
//    cout<<"++++++++++++pur\n"<<pur<<endl;
//
//
//    MatrixXd GurT = Gur.transpose();
//    Matrix3d H = (GurT*Gur).inverse();
////    cout<<"+++++H\n"<<H<<endl;
//
//    Vector3d bur = H*(GurT*pur);
//    cout<<"++++++bur\n"<<bur<<endl;
//
//    xyz = referBase+bur;
//    cout<<"++++++++++++basexyz\n"<<referBase<<endl;
//    cout<<"++++++++++++xyz\n"<<xyz<<endl;
//
//    gnss->isPositioned = true;
//
//    Solution pos(rcvtow,xyz,vxyz);
//    pos.Show("+++++++xyz-end");
//    gnss->AddPosRecord(pos);
//
//    fprintf(gnss->log,"xyz:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,xyz(0),xyz(1),xyz(2));
//    fprintf(gnss->log,"LLA:time = ,%.5f, pos = ,%.5f,%.5f,%.5f\n",rcvtow,LLA(1)*180/GPS_PI,LLA(0)*180/GPS_PI,LLA(2));
//
//    double pc = (svCectre->position - xyz).norm();
//    tui = svCectre->prMes+svCectre->tsDelta*Light_speed-svCectre->I-svCectre->T-pc;
//    gnss->tu = tui;
//    printf("tui= %f\n", tui);
//    char gga[128];
//    MakeGGA(gga,LLA,timeSolver);
//    rtk->SentGGA(gga,strlen(gga));
//    delete (this);
//}
int PosSolver::PositionRtk(vector<SV *> _svsIn) {
    int sigInd = 1;
    printf("rtk pr========= count=%d,nSat=%d,nSys=%d\n", gnss->count,nSat,nSys);
    PrepareSVsData(_svsIn);
    ProcessRtkData();

    if(nSat-nSys<4){
        printf("Not enough svs nSat,nSys=%d,%d\n",nSat,nSys);
        return -1;
    }
    Vector3d dtb,base = gnss->rtkManager.ECEF_XYZ;
    VectorXd pr_rb(nSat-nSys);
    MatrixXd G(nSat-nSys,3);
    int yhead = 0;

    for(SvSys* sys:svsBox.sysUsed){
        int Ni = sys->table.size();
        VectorXd rs(Ni);
        pr_rb.head(Ni-1)<<sys->DiffDouble(Diff_Pr_rb,sigInd);
        MatrixXd Ei = sys->GetE(base,rs);
        MatrixXd Di = sys->GetD();
        G.block(yhead,0,Ni-1,3)<<-Di*Ei;
        yhead+=Ni-1;
    }
    dtb = G.colPivHouseholderQr().solve(pr_rb);
    xyz = base+dtb;
    soltion = Solution(timeSolver,xyz,vxyz);
}


int PosSolverKalman::PositionKalman(vector<SV *> _svsIn) {
     int sigId = 1;
    printf("rtk pr========= count=%d,nSat=%d,nSys=%d\n", gnss->count,nSat,nSys);
    PrepareSVsData(_svsIn);
    ProcessRtkData();

    if(nSat-nSys<4){
        printf("Not enough svs nSat,nSys=%d,%d\n",nSat,nSys);
        return -1;
    }

    int ny = 2*(nSat-nSys),nx = nSat+6,yhead=0,xhead=6;
    Vector3d base = gnss->rtkManager.ECEF_XYZ;
    VectorXd x(nx),y(ny),hx(ny);
    MatrixXd P(nx),Ppred(nx),Q(nx);
    MatrixXd Hx(ny,nx),R(ny);

    x.fill(0);y.fill(0);hx.fill(0);P.fill(0);Q.fill(0);Hx.fill(0);R.fill(0);
    x.head(6)<<xyz,vxyz;
    P.block<6,6>(0,0)<<Pxv;
    Q.block<3,3>(0,0) = 1e5*Matrix3d::Identity();

    for(SvSys* sys:svsBox.sysUsed){
        double lambdai = GetFreq(sys->type,sigId,1);
        int Ni = sys->table.size();
        VectorXd rr(Ni),rb(Ni),Bi(Ni);
        Bi = sys->DiffZero(Diff_Cyc_0,sigId);
        sys->GetE(base,rb);

        MatrixXd Di = sys->GetD();
        MatrixXd DiT = Di.transpose();
        MatrixXd Ei = sys->GetE(xyz,rr);

        x.segment(xhead,Ni)<<Bi;

        y.segment(yhead,Ni-1)<<sys->DiffDouble(Diff_Cp_rb,1);
        y.segment(yhead+Ni-1,Ni-1)<<sys->DiffDouble(Diff_Pr_rb,1);

        VectorXd r_rb = Diff_Vec(rr-rb,Ni);
        VectorXd B_rb = Diff_Vec(Bi,Ni);
        hx.segment(yhead,Ni-1)<<r_rb+lambdai*B_rb;
        hx.segment(yhead+Ni-1,Ni-1)<<r_rb;

        Hx.block(yhead,0,Ni-1,3)<<-Di*Ei;
        Hx.block(yhead,xhead,Ni-1,Ni)<<lambdai*Di;
        Hx.block(yhead+Ni-1,0,Ni-1,3)<<-Di*Ei;

        VectorXd Pii = sys->DiffZero(Diff_cyP_0,sigId);
        P.block(xhead,xhead,Ni,Ni)=Pii.asDiagonal();

        VectorXd Rcpi = sys->DiffZero(Diff_Rcp_0,sigId);
        VectorXd Rpri = sys->DiffZero(Diff_Rpr_0,sigId);
        R.block(yhead,yhead,Ni-1,Ni-1)<<Di*Rcpi*DiT;
        R.block(yhead+Ni-1,yhead+Ni-1,Ni-1,Ni-1)<<Di*Rpri*DiT;

        yhead+=Ni-1;xhead+=Ni;
    }
}
int PosSolver::PositionRtkKalman() {
    Vector3d xyzBackup = xyz;
    int sigInd = 1;
    int result = 0;
    Vector3d base = gnss->rtkManager.ECEF_XYZ;
    vector<SV*> svsForCalcu[Nsys];
    PrepareSVsData(svsForCalcu);

    SolvePosition(svsForCalcu[7]);
    printf("ttttttttttttttttttttttttttttttttu %f\n", tui);

    ProcessRtkData(svsForCalcu);
    Vector3d xyzSingle = xyz;

    int N = 0,xi=6,yi=0;
    for (int sys = 0; sys < Nsys-1; ++sys)N+=sysEachNum[sys];
    printf("N= %d,nsys=%d\n", N,sysCount);
    if(N-sysCount<5){
        printf("svs not enough %d,nsys %d\n", N,sysCount);
        return -1;
    }
    VectorXd x(N+6),y(2*(N-sysCount)),hx(2*(N-sysCount));
    MatrixXd P(N+6,N+6),Q(N+6,N+6),Ppred(N+6,N+6);
    MatrixXd Hx(2*(N-sysCount),N+6),R(2*(N-sysCount),2*(N-sysCount));
    double cp_rb[2],pr_rb[2],r_rb[2],B_rb[2];
    vector<double*>pB,pPB,pRc,pRp;

    x.fill(0);y.fill(0);hx.fill(0);P.fill(0);Q.fill(0);Hx.fill(0);R.fill(0);
    x.head(6)<<xyzBackup,vxyz;
    P.block<6,6>(0,0)<<gnss->Pxv;
    Q.block<3,3>(0,0) = 1e5*Matrix3d::Identity();

    Vector3d bxlla,xyzrtk;
    //各系统
    for (int sys = 0,yHead,n,n_; sys < Nsys-1; ++sys) {
        n = sysEachNum[sys];
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
//                *tempBi = ((sv->cpMes*lambdai-sv->cpInterp[sigInd]) - (e_sr.norm()-e_sb.norm()) - tui)/lambdai;
//                sv->KalmanFirst = 0;
//            }

            double rdiff = sv->cpMes*lambdai-sv->prMes;
            double bdiff = sv->cpInterp[sigInd]-sv->prInterp[sigInd];
            double pr_ = sv->cpMes*lambdai-sv->cpInterp[sigInd];
            double pr__ = pr_ - tui;
            double r_ =  e_sr.norm()-e_sb.norm();
            double Biii = (rdiff-bdiff)/lambdai;
            printf("-----------cp- r-,tui,bi = %.2f,%.2f,%.2f,%.2f, %.2f\n",pr_,pr__,r_,tui,Biii);
            double rr = e_sr.norm(),rb = e_sb.norm();
            double ts = sv->tsDelta*Light_speed;
//            double tutemp = sv->prMes-rr+ts-sv->I-sv->T;
            double tub = sv->prInterp[sigInd]+ts-sv->I-sv->T-e_sb.norm();
            fprintf(gnss->logDebug,"%f,%f\n",tui,tub);
            double p_r_r = sv->prMes-rr-tui+ts;
            double p_r_b = sv->prInterp[sigInd] - rb - tub +ts;
//            fprintf(gnss->logDebug,"%.6f,\n",tui);
            if(sv->svId==0){
//                fprintf(gnss->logDebug,"%.6f,",sv->prInterp[sigInd]+sv->tsDelta*Light_speed-e_sb.norm());
//                fprintf(gnss->logDebug,"%.6f,",tui);
//                fprintf(gnss->logDebug,"%.6f,",sv->prMes+sv->tsDelta*Light_speed-e_sr.norm());
//
//                fprintf(gnss->logDebug,"%.6f,",Biii);
//                fprintf(gnss->logDebug,"%.6f,",sv->cpInterp[sigInd]);

                fprintf(gnss->logDebug,"%.6f,%.6f,",sv->I,sv->T);
                fprintf(gnss->logDebug,"%.6f,%.6f,",p_r_r,p_r_b);
                fprintf(gnss->logDebug,"%.6f,",tub);
                fprintf(gnss->logDebug,"%.6f,",tui);
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
        Solution pos00(rcvtow,xyzrtk,vxyz);
        pos00.Show("+++++xyz 111111 ");

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
        Solution pos(rcvtow,xyz,vxyz);
        pos.Show("+++++++xyz-222222");
        gnss->AddPosRecord(pos);
        cout<<"++++xyzDiff single\n"<<-gnss->xyz00<<endl;
        cout<<"++++xyzDiff rtk\n"<<xyzrtk-gnss->xyz00<<endl;
        cout<<"++++xyzDiff rtk kalman\n"<<xyzPred.head(3)-gnss->xyz00<<endl;

}

