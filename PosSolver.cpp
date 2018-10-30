#include "PosSolver.h"
#include "GNSS.h"
#include "EphemSp3.h"
#include "Svs.h"

auto Diff_Vec = [](VectorXd x,int n)->VectorXd { return x(0)*VectorXd::Ones(n-1)-x.tail(n-1);};

auto Diff_Pr_rb = [](SV* sv,int sigId)->double { return sv->measureDat.front()->prMes-sv->prInterp[sigId];};
auto Diff_Cp_rb = [](SV* sv,int sigId)->double { return sv->measureDat.front()->cpMes-sv->cpInterp[sigId];};
auto Diff_Cyc_0 = [](SV* sv,int sigId)->double& { return sv->measureDat.front()->cycle;};
auto Diff_cyP_0 = [](SV* sv,int sigId)->double& { return sv->measureDat.front()->cycleP;};
auto Diff_Rcp_0 = [](SV* sv,int sigId)->double { return 2*pow(sv->measureDat.front()->stdevCp,2);};
auto Diff_Rpr_0 = [](SV* sv,int sigId)->double { return 2*pow(sv->measureDat.front()->stdevPr,2);};


PosSolver::PosSolver(){}

PosSolver::PosSolver(SvAll svs, NtripRTK *rtk, GNSS *gnss):svsBox(svs),rtk(rtk),gnss(gnss){
    if(gnss->records.size()>0){
        Solution recent = gnss->records.front();
        xyz = recent.xyz;
        memcpy(tu,recent.tu,Nsys* sizeof(double));
    } else{
        xyz = gnss->xyzDefault;
    }
    XYZ2LLA(xyz,lla);
}

PosSolver::~PosSolver(){
    for(SvSys* sys:svsBox.sysUsed)delete(sys);
}

int PosSolver::PrepareSVsData(vector<SV*> &svsIn) {
    svsBox = gnss->svsManager;
    XYZ2LLA(xyz,lla);
    SelectSvsFromVisible(svsIn);
    if(svsBox.svUsedAll.empty()){ return -1;}
    timeSolver = svsBox.svUsedAll.front()->measureDat.front()->time;
    UpdateSvsPosition(svsBox.svUsedAll, timeSolver, gnss->ephemType);

}

int PosSolver::PositionSingle(vector<SV*> _svsIn) {
    PrepareSVsData(_svsIn);
    printf("single sinlel========== count=%d,nSat=%d,nSys=%d\n", gnss->count,nSat,nSys);

    soltion = Solution(timeSolver,xyz,vxyz,tu);
    soltion.Show("0000000before");
    if(nSat-nSys<4){
        printf("Not enough svs nSat,nSys=%d,%d\n",nSat,nSys);
        return -1;
    }
    double threshold = 0.1;
    VectorXd dt_xyz_tu(3+nSys),b(nSat);
    MatrixXd H(3+nSys,3+nSys);
    MatrixXd data(nSat,5);
    dt_xyz_tu.fill(1);
    int i =0,count =0;

    auto CalcuPc = [](SV* sv,int sigId)->double{
        return sv->measureDat.front()->prMes + Light_speed * sv->tsDelta - sv->I - sv->T;
    };
    auto CalcuI = [](SV* sv,int sigId)->double{ return  sv->I; };
    auto CalcuT = [](SV* sv,int sigId)->double{ return  sv->T; };
    auto Calcue = [](SV* sv,int sigId)->double{ return  sv->elevationAngle; };

    while (dt_xyz_tu.norm()>threshold){
        int yhead=0,xhead=3;
        MatrixXd G(nSat,3+nSys),Gt(3+nSys,nSat);
        G.fill(0);
        for(SvSys* sys:svsBox.sysUsed){
            int Ni = sys->table.size();
            if(0==Ni)continue;
            VectorXd rr(Ni),bi(Ni),pci;
            MatrixXd Ei = sys->GetE(xyz,rr);
            pci = sys->DiffZero(CalcuPc,1);
            bi = pci-rr-tu[sys->type]*(VectorXd::Ones(Ni));

            sys->MakeDebug(2);
            sys->AddAnaData(bi);
//            sys->Show();

            b.segment(yhead,Ni)<<bi;
            G.block(yhead,0,Ni,3)<<-Ei;
            G.block(yhead,xhead,Ni,1)<<VectorXd::Ones(Ni);
            yhead+=Ni;xhead++;
        }
//        dt_xyz_tu =  G.colPivHouseholderQr().solve(b);
        Gt = G.transpose();
        H = (Gt*G).inverse();
        dt_xyz_tu = H*Gt*b;
        xyz+=dt_xyz_tu.head(3);
//        cout<<"dt = "<<endl<<dt_xyz_tu<<endl;
//        cout<<"b = "<<endl<<b<<endl;
//        cout<<"G= "<<endl<<G<<endl;
//        cout<<"data = "<<endl<<data<<endl;

        i=3;
        //todo:sysused:
        for(SvSys* sys:svsBox.sysUsed)
            if(!sys->table.empty())
                tu[sys->type]+=dt_xyz_tu(i++);
        if(++count>10)return -1;
    }
    PDOP = sqrt(H(0,0)+H(1,1)+H(2,2));
    soltion = Solution(timeSolver,xyz,vxyz,tu);
    return 0;
}

int PosSolver::SelectSvsFromVisible(vector<SV*> &all) {
    for(SV* sv:all)
    {
        Measure *msTemp = sv->measureDat.front();
//        printf("getsv:%d,%02d:h:%d,pr=%lf,cp=%lf\t\t",sv->type,sv->svId,sv->SatH1,
//                msTemp->prMes,msTemp->cpMes);

        //1,is used?
        if(!sv->open){
            sprintf(sv->tip,"closed");
            continue;
        }
        //2,judge ephemeric
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
       sprintf(sv->tip,"++ok");
    }
    for(SV* sv:all) sv->PrintInfo(0);
    nSat = svsBox.svUsedAll.size();
    nSys = svsBox.sysUsed.size();
}

int PosSolver::UpdateSvsPosition(vector<SV *> &svs, GnssTime rt, int ephType) {
    printf("N0= %d\n", svs.size());
    printf("rcvtow %.3f\n", rt.tow);
    double ep[6];
    rt.time2epoch(ep);
    for(SV* sv:svs){
        Measure *measure = sv->measureDat.front();
        double timeForECEF = rt.tow;
        if(SYS_BDS == sv->type)timeForECEF-=14;
        timeForECEF-=tu[sv->type]/Light_speed;

        Sp3Cell cell;
        GnssTime ts0 = rt;
        switch (ephType){
            case 1:
                ts0+=(-measure->prMes/Light_speed-sv->tsDelta);
                if(EphemSp3::Sp32ECEF(sv->ephemSp3->records, ts0, cell))break;
                printf("ts2 =  %.3f\n", cell.ts*1e9);
                sv->xyz = cell.pxyz;
                break;
            case 0:
                sv->CalcuTime(timeForECEF);
                sv->CalcuECEF(sv->tsReal);
                break;
            default:
                continue;
        }

        double tt = ((sv->xyz-xyz).norm()+sv->I+sv->T)/Light_speed;
        EarthRotate(sv->xyz,sv->xyzR,tt);
        XYZ2LLA(sv->xyz,sv->lla);
        XYZ2LLA(sv->xyzR,sv->llaR);
        sv->CorrectIT(xyz,lla,timeForECEF);
    }
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
//    SvAll tempAll = gnss->svsManager;
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
//        sys->Show();
        int Ni = sys->table.size();
        VectorXd rr(Ni),pr_rbi(Ni-1);
        pr_rbi<<sys->DiffDouble(Diff_Pr_rb,sigInd);
        MatrixXd Ei = sys->GetE(base,rr);
        MatrixXd Di = sys->GetD();
        G.block(yhead,0,Ni-1,3)<<-Di*Ei;
        pr_rb.segment(yhead,Ni-1)<<pr_rbi;
        yhead+=Ni-1;

        VectorXd dtbi = (-Di*Ei).colPivHouseholderQr().solve(pr_rbi);
//        cout<<"Gi= "<<endl<<-Di*Ei<<endl;
//        cout<<"prrb_i= "<<endl<<pr_rbi<<endl;
//        cout<<"dtbi"<<endl<<dtbi<<endl;
//        cout<<"delta "<<endl<<pr_rbi-(-Di*Ei*dtbi);
        sys->MakeDebug(2);
//        sys->Show();
    }
//    cout<<"G= "<<endl<<G<<endl;
//    cout<<"prrb= "<<endl<<pr_rb<<endl;
    dtb = G.colPivHouseholderQr().solve(pr_rb);
    cout<<"dtb= "<<endl<<dtb<<endl;
    xyz = base+dtb;
    soltion = Solution(timeSolver,xyz,vxyz,tu);
}


int PosSolver::PositionKalman(vector<SV *> _svsIn) {
    int sigId = 1;
    svsBox = gnss->svsManager;
    PrepareSVsData(_svsIn);
    ProcessRtkData();

    printf("rtk pr========= count=%d,nSat=%d,nSys=%d\n", gnss->count,nSat,nSys);

    if(nSat-nSys<4){
        printf("Not enough svs nSat,nSys=%d,%d\n",nSat,nSys);
        return -1;
    }

    int ny = 2*(nSat-nSys),nx = nSat+6,yhead=0,xhead=6;
    Vector3d base = gnss->rtkManager.ECEF_XYZ;
    VectorXd x(nx),y(ny),hx(ny);
    MatrixXd P(nx,nx),Ppred(nx,nx),Q(nx,nx);
    MatrixXd Hx(ny,nx),R(ny,ny),G(nSat-nSys,nx);

    x.fill(0);y.fill(0);hx.fill(0);P.fill(0);Q.fill(0);Hx.fill(0);R.fill(0),G.fill(0);
    x.head(6)<<xyz,vxyz;
    P.block<6,6>(0,0)<<Pxv;
    Q = MatrixXd::Identity(nx,nx);
//    Q.block<3,3>(0,0) = 1e5*Matrix3d::Identity();
//    准备数据
    for(SvSys* sys:svsBox.sysUsed){
        double lambdai = GetFreq(sys->type,sigId,1);
        int Ni = sys->table.size();
        VectorXd rr(Ni),rb(Ni),Bi(Ni);
        Bi = sys->DiffZero(Diff_Cyc_0,sigId);
        sys->GetE(base,rb);

        MatrixXd Di = sys->GetD();
        MatrixXd DiT = Di.transpose();
        MatrixXd Ei = sys->GetE(xyz,rr);
        G.block(yhead/2,xhead,Ni-1,Ni)<<Di;

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
//        cout<<"Di"<<endl<<Di<<endl;
//        cout<<"Rcpi"<<endl<<Rcpi<<endl;
//        cout<<"R"<<endl<<R<<endl;
        R.block(yhead,yhead,Ni-1,Ni-1)<<Di*Rcpi.asDiagonal()*DiT;
        R.block(yhead+Ni-1,yhead+Ni-1,Ni-1,Ni-1)<<Di*Rpri.asDiagonal()*DiT;

        yhead+=2*(Ni-1);xhead+=Ni;


        sys->MakeDebug(2);
        sys->Show();
    }
    //计算Kalman
    P = P+Q;
    MatrixXd Ht = Hx.transpose();
    MatrixXd Kk = P*Ht*((Hx*P*Ht+R).inverse());
    Ppred = (MatrixXd::Identity(nx,nx)-Kk*Hx)*P;
    VectorXd xyzPred = x + Kk*(y-hx);

    //        cout<<"R"<<endl<<R<<endl;
    cout<<"H"<<endl<<Hx<<endl;

    cout<<"P"<<endl<<P<<endl;
    cout<<"HxPHt"<<endl<<P*Ht<<endl;
    cout<<"K"<<endl<<Kk<<endl;
    cout<<"Ppred"<<endl<<Ppred<<endl;
    cout<<"x"<<endl<<x.transpose()<<endl;
    cout<<(xyzPred-x).transpose()<<endl;
    cout<<xyzPred.transpose()<<endl;
    cout<<"B_double"<<endl<<(G*xyzPred).transpose()<<endl;
    cout<<"y"<<endl<<y.transpose()<<endl<<hx.transpose()<<endl;

    //更新参数
    xyz = xyzPred.head(3);
    vxyz = xyzPred.segment(3,3);
    Pxv = Ppred.block<6,6>(0,0);
    xhead=6;
    for(SvSys* sys:svsBox.sysUsed){
        int Ni = sys->table.size();
        VectorXd Bii = xyzPred.segment(xhead,Ni);
        VectorXd Pii = Ppred.diagonal().segment(xhead,Ni);
        sys->SetValue(Diff_Cyc_0,sigId,Bii);
        sys->SetValue(Diff_cyP_0,sigId,Pii);
        xhead+=Ni;
    }

    soltion = Solution(timeSolver,xyz,vxyz,tu);
    return 0;
}
int PosSolver::InitKalman(GNSS *_gnss) {
    gnss = _gnss;
    rtk = &(gnss->rtkManager);
    xyz = gnss->xyzDefault;
    Pxv = 1e9*MatrixXd::Identity(6,6);
}
int PosSolver::Init(GNSS *_gnss) {
    gnss = _gnss;
    rtk = &(gnss->rtkManager);
    xyz = gnss->xyzDefault;
}
//int PosSolver::PositionRtkKalman() {
//    Vector3d xyzBackup = xyz;
//    int sigInd = 1;
//    int result = 0;
//    Vector3d base = gnss->rtkManager.ECEF_XYZ;
//    vector<SV*> svsForCalcu[Nsys];
//    PrepareSVsData(svsForCalcu);
//
//    SolvePosition(svsForCalcu[7]);
//    printf("ttttttttttttttttttttttttttttttttu %f\n", tui);
//
//    ProcessRtkData(svsForCalcu);
//    Vector3d xyzSingle = xyz;
//
//    int N = 0,xi=6,yi=0;
//    for (int sys = 0; sys < Nsys-1; ++sys)N+=sysEachNum[sys];
//    printf("N= %d,nsys=%d\n", N,sysCount);
//    if(N-sysCount<5){
//        printf("svs not enough %d,nsys %d\n", N,sysCount);
//        return -1;
//    }
//    VectorXd x(N+6),y(2*(N-sysCount)),hx(2*(N-sysCount));
//    MatrixXd P(N+6,N+6),Q(N+6,N+6),Ppred(N+6,N+6);
//    MatrixXd Hx(2*(N-sysCount),N+6),R(2*(N-sysCount),2*(N-sysCount));
//    double cp_rb[2],pr_rb[2],r_rb[2],B_rb[2];
//    vector<double*>pB,pPB,pRc,pRp;
//
//    x.fill(0);y.fill(0);hx.fill(0);P.fill(0);Q.fill(0);Hx.fill(0);R.fill(0);
//    x.head(6)<<xyzBackup,vxyz;
//    P.block<6,6>(0,0)<<gnss->Pxv;
//    Q.block<3,3>(0,0) = 1e5*Matrix3d::Identity();
//
//    Vector3d bxlla,xyzrtk;
//    //各系统
//    for (int sys = 0,yHead,n,n_; sys < Nsys-1; ++sys) {
//        n = sysEachNum[sys];
//        n_=n-1;
//        if(0==n)continue;
//        cout <<"sys="<< sys <<", n="<<n<<endl;
//
//        double lambdai = GetFreq((SysType)sys, sigInd, true);
//        printf("lambdai %f\n", lambdai);
//        MatrixXd E(n,3),D(n-1,n),Rci(n,n),Rpi(n,n);
//        D.fill(0);Rci.fill(0);Rpi.fill(0);
//        //各卫星
//        VectorXd yp(n-1);
//
//        for (int i = 0,i_; i < n; ++i) {
//            SV *sv = svsForCalcu[sys][i];
//            Vector3d svPositionEarthRotate;
//            EarthRotate(sv->position,svPositionEarthRotate,(sv->position - xyz).norm()/Light_speed);
//            Vector3d e_sr = svPositionEarthRotate - xyz, e_sb = svPositionEarthRotate -base;
//            E.block<1,3>(i,0)<<(e_sr/e_sr.norm()).transpose();
//
//            printf("xi,yi %d,%d,\t", xi,yi);
//            printf("rtksv:%d,%02d:h:%d,pr=%lf---%lf,.,cp=%lf,%.1f---%lf,%.1f [%.2f]\telev=%lf\n",sv->type,sv->svId,sv->SatH1,
//                   sv->prMes,sv->prInterp[sigInd],sv->cpMes*lambdai,sv->cpMes-sv->prMes/lambdai,sv->cpInterp[sigInd],
//                   (sv->cpInterp[sigInd]-sv->prInterp[sigInd])/lambdai,sv->cpMes-sv->cpInterp[sigInd]/lambdai,sv->elevationAngle);
//            //单差部分
//            double *tempBi=&(gnss->cycle[sys][sv->svId]);     pB.push_back(tempBi);
//            double *tempPB=&(gnss->PB[sys][sv->svId]);        pPB.push_back(tempPB);
//            double *tempRci=&(gnss->sigmaCy[sys][sv->svId]);  pRc.push_back(tempRci);
//            double *tempRpi=&(gnss->sigmaPr[sys][sv->svId]);  pRp.push_back(tempRpi);
//
////            if(sv->KalmanFirst){
////                *tempBi = ((sv->cpMes*lambdai-sv->cpInterp[sigInd]) - (e_sr.norm()-e_sb.norm()) - tui)/lambdai;
////                sv->KalmanFirst = 0;
////            }
//
//            double rdiff = sv->cpMes*lambdai-sv->prMes;
//            double bdiff = sv->cpInterp[sigInd]-sv->prInterp[sigInd];
//            double pr_ = sv->cpMes*lambdai-sv->cpInterp[sigInd];
//            double pr__ = pr_ - tui;
//            double r_ =  e_sr.norm()-e_sb.norm();
//            double Biii = (rdiff-bdiff)/lambdai;
//            printf("-----------cp- r-,tui,bi = %.2f,%.2f,%.2f,%.2f, %.2f\n",pr_,pr__,r_,tui,Biii);
//            double rr = e_sr.norm(),rb = e_sb.norm();
//            double ts = sv->tsDelta*Light_speed;
////            double tutemp = sv->prMes-rr+ts-sv->I-sv->T;
//            double tub = sv->prInterp[sigInd]+ts-sv->I-sv->T-e_sb.norm();
//            fprintf(gnss->logDebug,"%f,%f\n",tui,tub);
//            double p_r_r = sv->prMes-rr-tui+ts;
//            double p_r_b = sv->prInterp[sigInd] - rb - tub +ts;
////            fprintf(gnss->logDebug,"%.6f,\n",tui);
//            if(sv->svId==0){
////                fprintf(gnss->logDebug,"%.6f,",sv->prInterp[sigInd]+sv->tsDelta*Light_speed-e_sb.norm());
////                fprintf(gnss->logDebug,"%.6f,",tui);
////                fprintf(gnss->logDebug,"%.6f,",sv->prMes+sv->tsDelta*Light_speed-e_sr.norm());
////
////                fprintf(gnss->logDebug,"%.6f,",Biii);
////                fprintf(gnss->logDebug,"%.6f,",sv->cpInterp[sigInd]);
//
//                fprintf(gnss->logDebug,"%.6f,%.6f,",sv->I,sv->T);
//                fprintf(gnss->logDebug,"%.6f,%.6f,",p_r_r,p_r_b);
//                fprintf(gnss->logDebug,"%.6f,",tub);
//                fprintf(gnss->logDebug,"%.6f,",tui);
//                fprintf(gnss->logDebug,"%.6f,",rdiff);
//                fprintf(gnss->logDebug,"%.6f,",bdiff);
//                fprintf(gnss->logDebug,"%.6f,",(rdiff-bdiff)/lambdai);
//                fprintf(gnss->logDebug,"\n");
//            }
//            //            *tempBi = Biii;
//            x(xi) = *tempBi;
//            P(xi,xi) = *tempPB;
//            Rci(i,i) = 2*pow(*tempRci,2);
//            Rpi(i,i) = 2*pow(*tempRpi,2);
//
//            xi++;
//            cp_rb[1] = sv->cpMes*lambdai-sv->cpInterp[sigInd];
//            pr_rb[1] = sv->prMes-sv->prInterp[sigInd];
//            B_rb[1] = *tempBi;
//            r_rb[1] = e_sr.norm()-e_sb.norm();
//
//            //双差部分
//            if(0==i){
//                cp_rb[0] = cp_rb[1];    pr_rb[0] = pr_rb[1];
//                B_rb[0] = B_rb[1];      r_rb[0] = r_rb[1];
//
//                yHead = 2*yi;
//                continue;
//            }
//            i_ = i-1;
//            D(i_,0)=1;D(i_,i)=-1;
//            y(yHead+i_)     = cp_rb[0] - cp_rb[1];
//            y(yHead+n_+i_)  = pr_rb[0] - pr_rb[1];
//            yp(i_) = y(yHead+n_+i_);
//
//            hx(yHead+i_)    = r_rb[0]-r_rb[1] + lambdai*(B_rb[0] - B_rb[1]);
//            hx(yHead+n_+i_) = r_rb[0]-r_rb[1];
//
//            printf("yi comp [%.2f,%.2f]\n", y(yHead+i_), y(yHead+n_+i_));
//            printf("hxi comp [%.2f,%.2f]\n", hx(yHead+i_), hx(yHead+n_+i_));
//            printf("B0,bi %f --- %f\n", B_rb[0],*tempBi);
//            yi++;
//
//            if(sv->svId==0){
//                double rrb12 = r_rb[0]-r_rb[1];
//                double prb12 = pr_rb[0]- pr_rb[1];
//                fprintf(gnss->logDebug,"%f,%f,%f,%f\n",prb12,rrb12,prb12-rrb12,sv->tsDelta*Light_speed);
//            }
//        }
//        Vector3d bx = (-D*E).colPivHouseholderQr().solve(yp);
//        xyzrtk = bx+base;
////        x.head(3) = xyzrtk;
//        Solution pos00(rcvtow,xyzrtk,vxyz);
//        pos00.Show("+++++xyz 111111 ");
//
////        cout<<D<<E<<endl;
//        MatrixXd Dt = D.transpose();
//        Hx.block(yHead,0,n-1,3)<<-D*E;
//        Hx.block(yHead+n-1,0,n-1,3)<<-D*E;
//        Hx.block(yHead,xi-n,n-1,n)<<lambdai*D;
//        R.block(yHead,yHead,n-1,n-1)<<D*Rci*Dt;
//        R.block(yHead+n-1,yHead+n-1,n-1,n-1)<<D*Rpi*Dt;
//    }
//    P = 1*MatrixXd::Identity(N+6,N+6)+P;
////    P.block<3,3>(0,0) = P.block<3,3>(0,0)+0.1*Matrix3d::Identity();
//
//    MatrixXd Ht = Hx.transpose();
//    MatrixXd Kk = P*Ht*((Hx*P*Ht+R).inverse());
//    Ppred = (MatrixXd::Identity(N+6,N+6)-Kk*Hx)*P;
//    VectorXd xyzPred = x + Kk*(y-hx);
//
//    cout <<"\nx\n"<<x.transpose() <<"\ny\n"<<y.transpose()<<"\nhx\n"<<hx.transpose()<<"\ny-hx\n"<<(y-hx).transpose()<<endl;
////    cout<<"\nP\n"<<P <<"\nHx\n"<<Hx<<"\npPred\n"<<Ppred<<endl;
////    cout <<"\nKk\n"<<Kk<<endl;
////    cout <<"\nR\n"<<R<<endl<<;
//    cout<<"\nxyzAdd\n"<<(Kk*(y-hx)).transpose() <<endl;
//    cout<<"\nxyzPred\n"<<xyzPred.transpose() <<endl;
//    //更新
//    xyz = xyzPred.head(3);
////    xyz = xyzrtk;
//    vxyz = xyzPred.segment(3,3);
//    for (int i = 0; i < pB.size(); ++i) {
//        *pB[i] = xyzPred(6+i);
//        *pPB[i] = Ppred(6+i,6+i);
//    }
//    gnss->Pxv = Ppred.block<6,6>(0,0);
//        //    gnss->Pxv(0,0) = Ppred(0,0);
//        //    gnss->Pxv(1,1) = Ppred(1,1);
//    //    gnss->Pxv(2,2) = Ppred(2,2);
//
//        //显示
//        gnss->isPositioned = true;
//        Solution pos(rcvtow,xyz,vxyz);
//        pos.Show("+++++++xyz-222222");
//        gnss->AddPosRecord(pos);
//        cout<<"++++xyzDiff single\n"<<-gnss->xyz00<<endl;
//        cout<<"++++xyzDiff rtk\n"<<xyzrtk-gnss->xyz00<<endl;
//        cout<<"++++xyzDiff rtk kalman\n"<<xyzPred.head(3)-gnss->xyz00<<endl;
//
//}

int PosSolver::AnaData(vector<SV *> _svsIn) {
    for(SV* sv:_svsIn){
        if(sv->type==SYS_GPS&&sv->svId==10){
            Measure* mes = sv->measureDat.front();
            fprintf(gnss->logDebug,"%.4f,%f,%f,%f,%f,",mes->time.tow,mes->prMes-mes->cpMes,mes->cpMes,mes->lockTime,mes->cno);
            fprintf(gnss->logDebug,"%f,%f,%f,%f,%f,%d\n",mes->stdevPr,mes->stdevCp,mes->stdevDo,mes->cycle,mes->cycleP,mes->trkStat);
        }
    }
}