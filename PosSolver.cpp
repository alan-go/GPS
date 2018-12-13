#include "PosSolver.h"
#include "GNSS.h"
#include "EphemSp3.h"
#include "EphemBst.h"
#include "Svs.h"

auto Diff_Vec = [](VectorXd x,int n)->VectorXd { return x(0)*VectorXd::Ones(n-1)-x.tail(n-1);};

auto Diff_Pr_rb = [](SV* sv,int sigId)->double { return sv->measureDat.front()->prMes-sv->prInterp[sigId];};
auto Diff_Cp_rb = [](SV* sv,int sigId)->double { return sv->measureDat.front()->cpMes-sv->cpInterp[sigId];};
auto Diff_Cyc_0 = [](SV* sv,int sigId)->double& {
    Measure* ms = sv->measureDat[0];
//    ms->cycle+=ms->cycleSlip;
    return ms->cycle;
};
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
    cout<<"---Prepare Data\n"<<endl;
    svsBox=SvAll();
    SelectSvsFromVisible(svsIn);
    if(svsBox.svUsedAll.empty()){
        cout<<"---Warning: Empty sv lists\n"<<endl;
        return -1;
    }
    timeSolLast = timeSol;
    timeSol = gnss->rTime;
    memcpy(tu,gnss->solSingle.tu,Nsys* sizeof(double));
    UpdateSvsPosition(svsBox.svUsedAll, timeSol, gnss->ephemType);
    nSat = svsBox.UpdateSVs("elev");
    nSys = svsBox.sysUsed.size();
    return 0;
}

int PosSolver::PositionSingle(vector<SV*> _svsIn) {
//    PrepareSVsData(_svsIn);
    printf("single sinlel========== count=%d,nSat=%d,nSys=%d\n", gnss->count,nSat,nSys);

    solSingle = Solution(timeSol,xyz,vxyz,tu);
    solSingle.Show("0000000before");
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
        return sv->measureDat.front()->prMes + Light_speed * sv->tsdt - sv->I - sv->T;
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
//
//            sys->MakeDebug(2);
//            sys->AddAnaData(bi);
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
    solSingle = Solution(timeSol,xyz,vxyz,tu);
    ///////////////////////DeBug/////////////////////////
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

            VectorXd rrRAC(Ni),biRAC(Ni),pciRAC;
            MatrixXd Eir = sys->GetE(gnss->xyzRAC,rrRAC);
            pciRAC = sys->DiffZero(CalcuPc,1);
            biRAC = pci-rrRAC-tu[sys->type]*(VectorXd::Ones(Ni));

            sys->MakeDebug(2);
            sys->AddAnaData(bi);
            sys->AddAnaData(biRAC);
            sys->Show();

            b.segment(yhead,Ni)<<bi;
            G.block(yhead,0,Ni,3)<<-Ei;
            G.block(yhead,xhead,Ni,1)<<VectorXd::Ones(Ni);
            yhead+=Ni;xhead++;
        }
    ///////////////////////DeBug/////////////////////////
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
        if(!sv->CheckEphemStates(timeSol,gnss->ephemType))continue;
        //3,measure
//        sv->SmoothKalman0();
//        if(sv->kal.state<30)continue;
        if(!sv->MeasureGood())continue;
        //4,elevtion angle
//        if(!sv->ElevGood())continue;

        sv->DetectCycleSlip();

        svsBox.AddToUsed(sv);
       sprintf(sv->tip,"++ok");
    }
    for(SV* sv:all) sv->PrintInfo(0);
}

int PosSolver::UpdateSvsPosition(vector<SV *> &svs, GnssTime rTime, int ephType) {
    printf("\nN0= %d\n", svs.size());
    printf("rcvtow %.3f\n", rTime.tow);
//    double ep[6];
//    rTime.time2epoch(ep);
    for(SV* sv:svs){
        sv->PrintInfo(1);
        Measure *ms = sv->measureDat.front();
        double tow,tot;
        GnssTime ts0 = rTime;
        ts0 += - ms->prMes/Light_speed;
        //sp3文件使用的是GPS时(虽然使用年月日表示的)
        //todo 相对论修正
        switch (ephType){
            case 1:
                if(sv->ephemSp3->CalcuTs(ts0))return -1;
                if(sv->ephemSp3->CalcuECEF(sv->ts))return -1;
                break;
            case 0:
                if(SYS_BDS == sv->type)ts0+=-14;
                if(sv->ephemBst->CalcuTs(ts0))return -1;
                if(sv->ephemBst->CalcuECEF(sv->ts))return -1;

//                sv->CalcuTs(ts0);
//                sv->CalcuECEF(sv->ts);
                //todo tsDelta 换成米单位
                break;
            default:
                break;
        }
        sv->PrintInfo(1);
        tot = (ms->prMes-tu[sv->type])/Light_speed+sv->tsdt;
        sv->RotateECEF(tot);
        sv->CorrectIT(xyz,lla,tow);
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
//    SvAll temp=svsBox;
//    svsBox.sysUsed.clear();

    for(auto ites=svsBox.sysUsed.begin();ites!=svsBox.sysUsed.end();){
        SvSys* sys = *ites;
//        for(SvSys* sys:temp.sysUsed){
        for(auto ite=sys->table.begin();ite!=sys->table.end();){
            SV* sv = *ite;
            double time = timeSol.tow - tu[sv->type]/Light_speed;
            if(sv->type==SYS_BDS)time-=14;
            double pr = sv->InterpRtkData(time,sigInd);
//            if(0!=pr)svsBox.AddToUsed(sv);
            if(0==pr)ite=sys->table.erase(ite); else ite++;
        }
        if(sys->table.empty())ites=svsBox.sysUsed.erase(ites); else ites++;
    }
    nSat = svsBox.UpdateSVs();
    nSys = svsBox.sysUsed.size();
}

int PosSolver::PositionRtk(vector<SV *> _svsIn) {
    int sigInd = 1;
    printf("rtk pr========= count=%d,nSat=%d,nSys=%d\n", gnss->count,nSat,nSys);
//    PrepareSVsData(_svsIn);
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
        VectorXd rr(Ni),pr_rbi(Ni-1);
        pr_rbi<<sys->DiffDouble(Diff_Pr_rb,sigInd);
        MatrixXd Ei = sys->GetE(base,rr);
        MatrixXd Di = sys->GetD();
        G.block(yhead,0,Ni-1,3)<<-Di*Ei;
        pr_rb.segment(yhead,Ni-1)<<pr_rbi;
        yhead+=Ni-1;

        VectorXd dtbi = (-Di*Ei).colPivHouseholderQr().solve(pr_rbi);
        ///////////////////////DEbug/////////////////////////////
//        cout<<"Gi= "<<endl<<-Di*Ei<<endl;
//        cout<<"prrb_i= "<<endl<<pr_rbi<<endl;
//        cout<<"dtbi"<<endl<<dtbi<<endl;
//        cout<<"delta "<<endl<<pr_rbi-(-Di*Ei*dtbi);
        sys->MakeDebug(2);
//        sys->AddAnaData(pr_rbi-)
        sys->Show();

        ///////////////////////DEbug/////////////////////////////
    }
//    cout<<"G= "<<endl<<G<<endl;
//    cout<<"prrb= "<<endl<<pr_rb<<endl;
    dtb = G.colPivHouseholderQr().solve(pr_rb);
    cout<<"dtb= "<<endl<<dtb<<endl;

//    cout<<"rtk_pr_rb= "<<endl<<pr_rb<<endl;
//    cout<<"rtk_G*dtb= "<<endl<<G*dtb<<endl;
    cout<<"rtk_res1= "<<endl<<pr_rb-G*dtb<<endl<<endl;
    cout<<"rtk_res2= "<<endl<<pr_rb-G*(gnss->xyzRAC-base)<<endl<<endl;
    xyz = base+dtb;
    solSingle = Solution(timeSol,xyz,vxyz,tu);
}

int PosSolver::PositionKalman(vector<SV *> _svsIn) {
    int sigInd = 1;
    double lastTime = gnss->solKalmans[0].time.tow;
    double dt = timeSol.tow - lastTime;
    ProcessRtkData();

    printf("Kalman Info========= count=%d,nSat=%d,nSys=%d,dt=%f\n", gnss->count, nSat, nSys, dt);
    ShowV3(xyz,"XYZ before Kal1");


    if (nSat - nSys < 4) {
        printf("Not enough svs nSat,nSys=%d,%d\n", nSat, nSys);
        return -1;
    }

    int ny = 2 * (nSat - nSys), nx = nSat + 6, yhead = 0, xhead = 6;
    Vector3d base = gnss->rtkManager.ECEF_XYZ;
    VectorXd x(nx), y(ny), hx(ny);
    MatrixXd P(nx, nx), Ppred(nx, nx), Q(nx, nx);
    MatrixXd Hx(ny, nx), R(ny, ny), G(nSat - nSys, nx);

    x.fill(0);    y.fill(0);    hx.fill(0);
    P.fill(0);    Q.fill(0);    Hx.fill(0);
    R.fill(0), G.fill(0);
    x.head(6) << xyz, vxyz;
    P.block<6, 6>(0, 0) << P66;
    Q = MatrixXd::Identity(nx, nx);
    Q.block(0, 0, 3, 3) *= pow(10 * dt, 2);
    Q.block(6, 6, nSat, nSat) *= 0.5;
//    Q.block<3,3>(0,0) = 1e5*Matrix3d::Identity();
//    准备数据
    for (SvSys *sys:svsBox.sysUsed) {
        double lambdai = GetFreq(sys->type, sigInd, 1);
        int Ni = sys->table.size();
        VectorXd rr(Ni), rb(Ni), Bi(Ni);
//        Bi = sys->DiffZero(Diff_Cyc_0,sigInd);
        for (int i = 0; i < Ni; i++) {
            Measure *ms0 = sys->table[i]->measureDat[0];
            Measure *ms1 = sys->table[i]->measureDat[1];
            double dt01 = ms0->time.tow - ms1->time.tow;
            ms0->cycle += ms0->cycleSlip;
            ms0->cycleP+=ms0->cycleSlipQ;
            if (ms0->trackTime < 1) {
//                ms0->cycle=0;
//                ms0->cycleP = pow(500*dt01*dt01, 2);
                ms0->cycleP += pow(10*dt01*dt01, 2);
                if (abs(ms0->cycleRes) > 0.10){
                    printf("cycleQ= %f\n", ms0->cycleRes);
                    ms0->cycleP += pow(1000*ms0->stdevDo*ms0->cycleRes,2);
                }
            }
            Bi(i) = ms0->cycle;
        }

        sys->GetE(base, rb);

        MatrixXd Di = sys->GetD();
        MatrixXd DiT = Di.transpose();
        MatrixXd Ei = sys->GetE(xyz, rr);
        G.block(yhead / 2, xhead, Ni - 1, Ni) << Di;

        x.segment(xhead, Ni) << Bi;

        y.segment(yhead, Ni - 1) << sys->DiffDouble(Diff_Cp_rb, 1);
        y.segment(yhead + Ni - 1, Ni - 1) << sys->DiffDouble(Diff_Pr_rb, 1);

        VectorXd r_rb = Diff_Vec(rr - rb, Ni);
        VectorXd B_rb = Diff_Vec(Bi, Ni);
        hx.segment(yhead, Ni - 1) << r_rb + lambdai * B_rb;
        hx.segment(yhead + Ni - 1, Ni - 1) << r_rb;

        Hx.block(yhead, 0, Ni - 1, 3) << -Di * Ei;
        Hx.block(yhead, xhead, Ni - 1, Ni) << lambdai * Di;
        Hx.block(yhead + Ni - 1, 0, Ni - 1, 3) << -Di * Ei;

        VectorXd Pii = sys->DiffZero(Diff_cyP_0, sigInd);
        P.block(xhead, xhead, Ni, Ni) = Pii.asDiagonal();

        VectorXd Rcpi = sys->DiffZero(Diff_Rcp_0, sigInd);
        VectorXd Rpri = sys->DiffZero(Diff_Rpr_0, sigInd);
//        cout<<"Di"<<endl<<Di<<endl;
//        cout<<"Rcpi"<<endl<<Rcpi<<endl;
//        cout<<"R"<<endl<<R<<endl;
        R.block(yhead, yhead, Ni - 1, Ni - 1) << Di * Rcpi.asDiagonal() * DiT;
        R.block(yhead + Ni - 1, yhead + Ni - 1, Ni - 1, Ni - 1) << Di * Rpri.asDiagonal() * DiT;

        yhead += 2 * (Ni - 1);
        xhead += Ni;
        ///////////////////////DEbug/////////////////////////////
//        for(SV* sv:sys->table){
//            Measure *ms = sv->measureDat.front();
//            double p_f_r = ms->cpMes-ms->prMes;
//            double p_f_b = sv->cpInterp[sigInd]-sv->prInterp[sigInd];
//            sv->FPrintInfo(0);
//            fprintf(sv->fpLog,"%f,%f,%f,%f\n",timeSol.tow,p_f_r,p_f_b,ms->cycle);
//        }
//        sys->MakeDebug(2);
//        sys->Show();
        ///////////////////////DEbug/////////////////////////////
    }
    //计算Kalman
    P = P + Q;
    MatrixXd Ht = Hx.transpose();
    MatrixXd Kk = P * Ht * ((Hx * P * Ht + R).inverse());
    Ppred = (MatrixXd::Identity(nx, nx) - Kk * Hx) * P;
    VectorXd xyzPred = x + Kk * (y - hx);

    cout<<"R"<<endl<<R<<endl;
    cout << "H" << endl << Hx << endl;

    cout << "P" << endl << P << endl;
//    cout << "HxPHt" << endl << P * Ht << endl;
//    cout << "K" << endl << Kk << endl;
//    cout << "Ppred" << endl << Ppred << endl;
    cout << "x" << endl << x.transpose() << endl;
    cout << (xyzPred - x).transpose() << endl;
    cout << xyzPred.transpose() << endl;
    cout << "B_double" << endl << (G * xyzPred).transpose() << endl;
    cout << "y" << endl << y.transpose() << endl << hx.transpose() << endl<< (y-hx).transpose() << endl;

    //更新参数
    xyz = xyzPred.head(3);
    vxyz = xyzPred.segment(3, 3);
    P66 = Ppred.block<6, 6>(0, 0);
    xhead = 6;
    for (SvSys *sys:svsBox.sysUsed) {
        int Ni = sys->table.size();
        VectorXd Bii = xyzPred.segment(xhead, Ni);
        VectorXd Pii = Ppred.diagonal().segment(xhead, Ni);
        sys->SetValue(Diff_Cyc_0, sigInd, Bii);
        sys->SetValue(Diff_cyP_0, sigInd, Pii);
        xhead += Ni;
    }

    solSingle = Solution(timeSol, xyz, vxyz, tu);
    return 0;
}

int PosSolver::ResetKalRtk(int N, int M, Kalman & kal, int L){
    static int n=0;
    cout<<"ResetKalSingle  "<<n++<<endl;
    kal.Alloc(N, M, L);
    if (kal.state == 0) {
        dts=0;
//        solKalDopp=gnss->solver.solSingle;
        solKalDopp=Solution(timeSol,gnss->xyzDefault,Vector3d(0,0,0),tu);
        memcpy(solKalDopp.tu,solSingle.tu,Nsys* sizeof(double));

        P99.block(0,0,3,3)=1e6*Matrix3d::Identity();
        P99.block(3,3,3,3)=1e3*Matrix3d::Identity();
        P99(6,6)=P99(7,7)=500;
        P99(8,8)=1e4;

        kal.AnnAdd99.block(0,3,3,3)=Matrix3d::Identity();
        kal.AnnAdd99(6,8)=1;
        kal.AnnAdd99(7,8)=1;
        kal.state=1;
    }
    for(SV* sv:svsBox.svUsedAll)if(sv->kalState==0){
        sv->Ii=sv->I;
        sv->IiP=pow((cos(sv->elevationAngle)+0.3)*30,2);
        Measure* ms = sv->measureDat[0];
        double lambda=GetFreq(sv->type,1,1);
//        ms->cycle=(ms->cpMes-ms->prMes+2*sv->I)/lambda;
        ms->cycleP=10000;
        sv->kalState=1;
    }

    kal.Pnn.block(0,0,9,9)=P99;
    kal.x.head(3)=xyz=solKalDopp.xyz;
    kal.Qnn.block(0,0,3,3)=1.12*Matrix3d::Identity();

    kal.x.segment(3,3)=vxyz=solKalDopp.vxyz;
    kal.Qnn.block(3,3,3,3)=0.31*Matrix3d::Identity();

    kal.x(6)=solKalDopp.tu[SYS_BDS];
    kal.x(7)=solKalDopp.tu[SYS_GPS];
    kal.x(8)=tuf;
    kal.Qnn.block(6,6,3,3)=0.21*Matrix3d::Identity();

}

int PosSolver::PositionKalman2(vector<SV *> _svsIn) {
    int sigInd = 1;
    double dts = timeSol.tow-timeSolLast.tow;
    ProcessRtkData();


    printf("Kalman Info========= count=%d,nSat=%d,nSys=%d\n", gnss->count,nSat,nSys);

    if(nSat-nSys<4){
        printf("Not enough svs nSat,nSys=%d,%d\n",nSat,nSys);
        return -1;
    }

    int ny = 3*nSat-2*nSys,nx = nSat+9,yhead,xhead;
    ResetKalRtk(nx,ny,kalRtk);
//    ShowV3(xyz,"XYZ before Kal2");

    Vector3d base = gnss->rtkManager.ECEF_XYZ;
//    读数据，x,
    yhead=0,xhead=9;
    for(SvSys* sys:svsBox.sysUsed){
        double lambdai = GetFreq(sys->type,sigInd,1);
        int Ni = sys->table.size();
        for(int i=0;i<Ni;i++){
            SV *sv = sys->table[i];
            Measure *ms0 = sv->measureDat[0];
            Measure *ms1 = sv->measureDat[1];
            double dt01 = ms0->time.tow - ms1->time.tow;
            ms0->cycle += ms0->cycleSlip;
            ms0->cycleP+=ms0->cycleSlipQ;
            if (ms0->trackTime < 1) {
                ms0->cycleP += pow(10 * dt01 * dt01, 2);
                if (abs(ms0->cycleRes) > 0.1){
                    ms0->cycleP += pow(1000*ms0->stdevDo*ms0->cycleRes,2);
                }
            }

            kalRtk.x(xhead + i) = ms0->cycle;
            kalRtk.Pnn(xhead + i, xhead + i) = ms0->cycleP;
//            kalRtk.Qnn(xhead + i, xhead + i) = ms0->cycleSlipQ;
        }
        xhead+=Ni;
    }
//    预测x,
//        cout<<"Pnnbefore:  "<<endl<<kalRtk.Pnn<<endl;
//        cout<<"Qnnbefore:  "<<endl<<kalRtk.Qnn<<endl;
   kalRtk.Predict(dts);
//    kalRtk.Pnn+=kalRtk.Qnn;
    xyz = kalRtk.x.head(3);
//    ShowV3(xyz,"Predict kal1");

    yhead=0,xhead=9;
//    读数据，y,
    for(SvSys* sys:svsBox.sysUsed) {
        double lambdai = GetFreq(sys->type, sigInd, 1);
        int Ni = sys->table.size();
        VectorXd rr(Ni), rb(Ni), Bi(Ni);
        Bi=kalRtk.x.segment(xhead,Ni);

        sys->GetE(base,rb);

        MatrixXd Di = sys->GetD();
        MatrixXd DiT = Di.transpose();
        MatrixXd Ei = sys->GetE(xyz,rr);

        VectorXd r_rb = Diff_Vec(rr-rb,Ni);
        VectorXd B_rb = Diff_Vec(Bi,Ni);

        kalRtk.y.segment(yhead,Ni-1)<<sys->DiffDouble(Diff_Cp_rb,1);
        kalRtk.hx.segment(yhead,Ni-1)<<r_rb+lambdai*B_rb;
        kalRtk.Hmn.block(yhead,0,Ni-1,3)<<-Di*Ei;
        kalRtk.Hmn.block(yhead,xhead,Ni-1,Ni)<<lambdai*Di;
        VectorXd Rcpi = sys->DiffZero(Diff_Rcp_0,sigInd);
        kalRtk.Rmm.block(yhead,yhead,Ni-1,Ni-1)<<Di*Rcpi.asDiagonal()*DiT;
        yhead+=Ni-1;

        kalRtk.y.segment(yhead,Ni-1)<<sys->DiffDouble(Diff_Pr_rb,1);
        kalRtk.hx.segment(yhead,Ni-1)<<r_rb;
        kalRtk.Hmn.block(yhead,0,Ni-1,3)<<-Di*Ei;
        VectorXd Rpri = sys->DiffZero(Diff_Rpr_0,sigInd);
        kalRtk.Rmm.block(yhead,yhead,Ni-1,Ni-1)<<Di*Rpri.asDiagonal()*DiT;
        yhead+=Ni-1;


        for (int i = 0; i < Ni; ++i) {
            SV *sv = sys->table[i];
            Vector3d ei = sv->xyzR - kalRtk.x.head(3);
            double ri = ei.norm();
            ei /= ri;
            Measure *ms0 = sv->measureDat[0];
            kalRtk.y(yhead+i)=ms0->doMes;
            kalRtk.hx(yhead+i)=ei.transpose()*(sv->vxyzR-vxyz)+tuf-sv->ephemBst->clk.a1*Light_speed;
            kalRtk.Rmm(yhead+i,yhead+i)=pow(ms0->stdevDo,2);
            kalRtk.Ri(yhead+i)=ms0->stdevDo;
            kalRtk.Hmn.block(yhead+i,3,1,3)=-ei.transpose();
            kalRtk.Hmn(yhead+i,8)=1;
        }
        yhead+=Ni;xhead+=Ni;

        ///////////////////////DEbug/////////////////////////////

        ///////////////////////DEbug/////////////////////////////
    }
    //校正
    kalRtk.Rectify();
    //更新参数
    xyz=kalRtk.x.head(3);
    vxyz=kalRtk.x.segment(3,3);
    tu[SYS_BDS]=kalRtk.x(6);
    tu[SYS_GPS]=kalRtk.x(7);
    tuf=kalRtk.x(8);
    P99=kalRtk.Pnn.block(0,0,9,9);

    xhead=9;
    for(SvSys* sys:svsBox.sysUsed) {
        int Ni = sys->table.size();
        for (int i = 0; i < Ni; i++) {
            SV *sv = sys->table[i];
            Measure *ms0 = sv->measureDat[0];
            ms0->cycle=kalRtk.x(xhead + i);
            ms0->cycleP=kalRtk.Pnn(xhead + i,xhead + i);
        }
        xhead+=Ni;
    }

    solKalDopp=Solution(timeSol,xyz,vxyz,tu);
    cout<<"KalRtkDopp"<<kalRtk.x.transpose()<<endl;
    return 0;
}


int PosSolver::InitKalman(GNSS *_gnss) {
    gnss = _gnss;
    rtk = &(gnss->rtkManager);
//    xyz = gnss->solSingle.xyz;
    xyz = gnss->xyzDefault;
    P66 = 1e6*MatrixXd::Identity(6,6);
}
int PosSolver::Init(GNSS *_gnss) {
    gnss = _gnss;
    rtk = &(gnss->rtkManager);
    xyz = gnss->xyzDefault;
}

int PosSolver::AnaData(vector<SV *> _svsIn) {
    for(SV* sv:_svsIn){
        if(sv->type==SYS_GPS&&sv->svId==10){
            Measure* mes = sv->measureDat.front();
            fprintf(gnss->logDebug,"%.4f,%f,%f,%f,%f,",mes->time.tow,mes->prMes-mes->cpMes,mes->cpMes,mes->lockTime,mes->cno);
            fprintf(gnss->logDebug,"%f,%f,%f,%f,%f,%d\n",mes->stdevPr,mes->stdevCp,mes->stdevDo,mes->cycle,mes->cycleP,mes->trkStat);
        }
    }
}

int PosSolver::PositionSingleNew(vector<SV *> _svsIn) {
    printf("single sinlel========== count=%d,nSat=%d,nSys=%d\n", gnss->count,nSat,nSys);

    solSingle = Solution(timeSol,xyz,vxyz,tu);
    solSingle.Show("0000000before");
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
        return sv->measureDat.front()->prCor + Light_speed * sv->tsdt- sv->I - sv->T;
    };
//    auto CalcuWi = [](SV* sv,int sigId)->double{ return  1/pow(sv->measureDat[0]->stdPrCor,2); };
    auto CalcuWi = [](SV* sv,int sigId)->double{ return  1; };
    auto CalcuI = [](SV* sv,int sigId)->double{ return  sv->I; };
    auto CalcuT = [](SV* sv,int sigId)->double{ return  sv->T; };
    auto Calcue = [](SV* sv,int sigId)->double{ return  sv->elevationAngle; };

    while (dt_xyz_tu.norm()>threshold){
        int yhead=0,xhead=3;
        MatrixXd G(nSat,3+nSys),Gt(3+nSys,nSat),W(nSat,nSat);
        G.fill(0);
        for(SvSys* sys:svsBox.sysUsed){
            int Ni = sys->table.size();
            if(0==Ni)continue;
            VectorXd rr(Ni),bi(Ni),wi,pci;
            MatrixXd Ei = sys->GetE(xyz,rr);
            pci = sys->DiffZero(CalcuPc,1);
            wi = sys->DiffZero(CalcuWi,1);
            bi = pci-rr-tu[sys->type]*(VectorXd::Ones(Ni));
//
            sys->MakeDebug(2);
            sys->AddAnaData(bi);
            sys->AddAnaData(wi);
            sys->Show();

            b.segment(yhead,Ni)<<bi;
            G.block(yhead,0,Ni,3)<<-Ei;
            G.block(yhead,xhead,Ni,1)<<VectorXd::Ones(Ni);
            W.block(yhead,yhead,Ni,Ni)=wi.asDiagonal();
            yhead+=Ni;xhead++;
        }
//        dt_xyz_tu =  G.colPivHouseholderQr().solve(b);
        MatrixXd WG = W*G;
        VectorXd Wb = W*b;
        MatrixXd WGT = WG.transpose();
        H = (WGT*WG).inverse();
        dt_xyz_tu = H*WGT*Wb;
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
    solSingle = Solution(timeSol,xyz,vxyz,tu);
    return 0;
}

int PosSolver::ResetKalSingle(int N, int M, vector<SV *> &svsIn, int L) {
    static int n=0;
    cout<<"ResetKalSingle  "<<n++<<endl;
    kalSingle.Alloc(N, M, L);
    if (kalSingle.state == 0) {
        dts=0;
        solKalSigle=solSingle;
        memcpy(solKalSigle.tu,solSingle.tu,Nsys* sizeof(double));

        P99.block(0,0,3,3)=1e6*Matrix3d::Identity();
        P99.block(3,3,3,3)=1e3*Matrix3d::Identity();
        P99(6,6)=P99(7,7)=500;
        P99(8,8)=1e4;

        kalSingle.AnnAdd99.block(0,3,3,3)=Matrix3d::Identity();
        kalSingle.AnnAdd99(6,8)=1;
        kalSingle.AnnAdd99(7,8)=1;
        kalSingle.state=1;
    }
    for(SV* sv:svsIn)if(sv->kalState==0){
        sv->Ii=sv->I;
        sv->IiP=pow((cos(sv->elevationAngle)+0.3)*30,2);
        Measure* ms = sv->measureDat[0];
        double lambda=GetFreq(sv->type,1,1);
        ms->cycle=(ms->cpMes-ms->prMes+2*sv->I)/lambda;
        ms->cycleP=5000;
        sv->kalState=1;
    }

    kalSingle.Pnn.block(0,0,9,9)=P99;
    kalSingle.x.head(3)=solKalSigle.xyz;
    kalSingle.Qnn.block(0,0,3,3)=0.2*Matrix3d::Identity();

    kalSingle.x.segment(3,3)=solKalSigle.vxyz;
    kalSingle.Qnn.block(3,3,3,3)=0.1*Matrix3d::Identity();

    kalSingle.x(6)=solKalSigle.tu[SYS_BDS];
    kalSingle.x(7)=solKalSigle.tu[SYS_GPS];
    kalSingle.x(8)=tuf;
    kalSingle.Qnn.block(6,6,3,3)=0.1*Matrix3d::Identity();

}

int PosSolver::PosKalSng(vector<SV *> _svsIn) {
    printf("\n\nResetKalSingle========== count=%d,nSat=%d,nSys=%d\n", gnss->count,nSat,nSys);
    int sigInd = 1;
    dts = timeSol.tow-timeSolLast.tow;
    int ny = 3*nSat,nx = 2*nSat+9,yhead=0,xhead=9;
    ResetKalSingle(nx,ny,_svsIn);

    yhead=0;xhead=9;
    for(SvSys* sys:svsBox.sysUsed) {
        double lambdai = GetFreq(sys->type, sigInd, 1);
        int Ni = sys->table.size();
        for (int i = 0; i < Ni; i++) {
            SV *sv = sys->table[i];
            Measure *ms0 = sv->measureDat[0];
            Measure *ms1 = sv->measureDat[1];
            double dt01 = ms0->time.tow - ms1->time.tow;
            ms0->cycle += ms0->cycleSlip;

            kalSingle.x(xhead + i) = ms0->cycle;
            kalSingle.Pnn(xhead + i, xhead + i) = ms0->cycleP;
            kalSingle.Qnn(xhead + i, xhead + i) = ms0->cycleSlipQ;
            kalSingle.x(xhead + Ni + i) = sv->Ii;
            kalSingle.Pnn(xhead + i + Ni, xhead + i + Ni) = sv->IiP;
            kalSingle.Qnn(xhead + i + Ni, xhead + i + Ni) = 0.5;
        }
        xhead+=2*Ni;yhead+=3*Ni;
    }

    kalSingle.Predict(dts);
    yhead=0;xhead=9;
     for(SvSys* sys:svsBox.sysUsed){
        double lambdai = GetFreq(sys->type,sigInd,1);
        int Ni = sys->table.size();
        for (int i = 0; i < Ni; ++i) {
            SV* sv = sys->table[i];
            int tui=7;
            if(SYS_BDS==sv->type)tui=6;
            Vector3d ei=sv->xyzR-kalSingle.x.head(3);
            double ri=ei.norm();
            ei/=ri;
            Measure *ms0 = sv->measureDat[0];
            Measure *ms1 = sv->measureDat[1];
            double dt01 = ms0->time.tow - ms1->time.tow;
            kalSingle.y(yhead+i)=ms0->cpMes;
            kalSingle.hx(yhead+i)=ri+kalSingle.x(tui)-sv->tsdt*Light_speed-sv->Ii+sv->T+lambdai*ms0->cycle;
            kalSingle.Rmm(yhead+i,yhead+i)=pow(ms0->stdevCp,2);
            kalSingle.Ri(yhead+i)=ms0->stdevCp;
            kalSingle.Hmn.block(yhead+i,0,1,3)=-ei.transpose();
            kalSingle.Hmn(yhead+i,tui)=1;
            kalSingle.Hmn(yhead+i,xhead+i)=lambdai;
            kalSingle.Hmn(yhead+i,xhead+i+Ni)=-1;

            kalSingle.y(yhead+i+Ni)=ms0->cpMes-ms0->prMes;
            kalSingle.hx(yhead+i+Ni)=-2*sv->Ii+lambdai*ms0->cycle;
            kalSingle.Rmm(yhead+i+Ni,yhead+i+Ni)=pow(ms0->stdevPr,2);
            kalSingle.Ri(yhead+i+Ni)=ms0->stdevPr;
            kalSingle.Hmn(yhead+i+Ni,xhead+i)=lambdai;
            kalSingle.Hmn(yhead+i+Ni,xhead+i+Ni)=-2;

            kalSingle.y(yhead+i+2*Ni)=ms0->doMes;
            kalSingle.hx(yhead+i+2*Ni)=ei.transpose()*(sv->vxyzR-vxyz)+tuf-sv->ephemBst->clk.a1*Light_speed;
            kalSingle.Rmm(yhead+i+2*Ni,yhead+i+2*Ni)=pow(ms0->stdevDo,2);
            kalSingle.Ri(yhead+i+2*Ni)=ms0->stdevDo;
            kalSingle.Hmn.block(yhead+i+2*Ni,3,1,3)=-ei.transpose();
            kalSingle.Hmn(yhead+i+2*Ni,8)=1;
        }
        xhead+=2*Ni;yhead+=3*Ni;
    }

    kalSingle.Rectify();
    xyz=kalSingle.x.head(3);
    vxyz=kalSingle.x.segment(3,3);
    tu[SYS_BDS]=kalSingle.x(6);
    tu[SYS_GPS]=kalSingle.x(7);
    tuf=kalSingle.x(8);
    P99=kalSingle.Pnn.block(0,0,9,9);

    yhead=0;xhead=9;
    for(SvSys* sys:svsBox.sysUsed) {
        int Ni = sys->table.size();
        for (int i = 0; i < Ni; i++) {
            SV *sv = sys->table[i];
            Measure *ms0 = sv->measureDat[0];
            ms0->cycle=kalSingle.x(xhead + i);
            ms0->cycleP=kalSingle.Pnn(xhead + i,xhead + i);
            sv->Ii=kalSingle.x(xhead + Ni + i);
            sv->IiP=kalSingle.Pnn(xhead + Ni + i,xhead + Ni + i);
        }
        xhead+=2*Ni;yhead+=3*Ni;
    }

    solKalSigle=Solution(timeSol,xyz,vxyz,tu);
    cout<<"KalSigX  "<<kalSingle.x.transpose()<<endl;
    return 0;
}