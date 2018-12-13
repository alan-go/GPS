#ifndef SVS_H
#define SVS_H

#include "CommonInclude.h"

using namespace std;
using namespace Eigen;


class GNSS;
class MSM4data;
class EphemSp3;
class EphemBst;
class SvAll;
struct Ionosphere{
    double a0,a1,a2,a3;
    double b0,b1,b2,b3;
    Ionosphere():a0(0),a1(0),a2(0),a3(0),b0(0),b1(0),b2(0),b3(0){}
};
class Kalman{
public:
    int state{0};//0:uninitialized
    int N,M,L;
    VectorXd x,y,hx,Ri;
    MatrixXd Pnn,Hmn,Rmm;
    MatrixXd Ann,AnnAdd99,Qnn;
    Kalman(){
        AnnAdd99=MatrixXd(9,9);AnnAdd99.fill(0);
    };
    Kalman(int n,int m,int l){
        Alloc(n,m,l);
//        Qnn(0,0)=50000;
    };
    int Alloc(int n,int m,int l){
        N=n;M=m;L=l;
        x=VectorXd(N);x.fill(0);
        y=VectorXd(M);y.fill(0);
        Ri=VectorXd(M);Ri.fill(0);
        hx=VectorXd(M);hx.fill(0);
        Pnn=MatrixXd(N,N);Pnn.fill(0);
        Hmn=MatrixXd(M,N);Hmn.fill(0);
        Rmm=MatrixXd(M,M);Rmm.fill(0);

        Ann=MatrixXd::Identity(N,N);
        Qnn=MatrixXd::Identity(N,N);Qnn.fill(0);
    }
    int Init(){
        Pnn = 10000*MatrixXd::Identity(N,N);
    };
    int Predict(double t){
        cout<<"x:  "<<endl<<x.transpose()<<endl;
        Ann.block(0,0,9,9)+=t*AnnAdd99;
//        cout<<"AnnAdd:  "<<endl<<Ann*x-x<<endl;
        x=Ann*x;
        Pnn=Ann*Pnn*Ann.transpose()+Qnn;

        cout<<"xPredict:  "<<endl<<x.transpose()<<endl;
//        cout<<"Ann:  "<<endl<<Ann<<endl;
//        cout<<"Qnn:  "<<endl<<Qnn<<endl;
//        cout<<"Pnn:  "<<endl<<Pnn<<endl;
    }
    int Rectify(){
        static int n=0;
        cout<<"Pnn:  "<<endl<<Pnn<<endl;
//        cout<<"Hmn:  "<<endl<<Hmn<<endl;
//        cout<<"Rmm:  "<<endl<<Rmm<<endl;
        MatrixXd Ht=Hmn.transpose();
        MatrixXd Kk=Pnn*Ht*(Hmn*Pnn*Ht+Rmm).inverse();

//        cout<<"Kk00:  "<<endl<<Hmn*Pnn*Ht+Rmm<<endl;
//        cout<<"Kk:  "<<endl<<Kk<<endl;


        cout<<"   x:  "<<x.transpose()<<endl;
        VectorXd xadd = Kk*(y-hx);
        x=x+xadd;
        Pnn = (MatrixXd::Identity(N,N)-Kk*Hmn)*Pnn;

        cout<<"xadd:  "<<xadd.transpose()<<endl;
        cout<<"   x:  "<<x.transpose()<<endl;
//        cout<<"  Ri:  "<<Ri.transpose()<<endl;
        cout<<"   y:  "<<y.transpose()<<endl;
        cout<<"  hx:  "<<(hx).transpose()<<endl;
        cout<<"y-hx:  "<<(y-hx).transpose()<<endl;
    };
};
class SV{
public:

    EphemSp3 *ephemSp3;
    EphemBst *ephemBst;
    int trackCount;
    bool ephemOkBst{0},ephemOkSp3{0};
    bool isBeiDouGEO;
    bool open{1},healthy{1}, measureGood{0}, elevGood{0};
    SysType type;
    int svId;///1,2,3...

    double I,T;
    Kalman kal;
    int kalState{0};
    FILE* fpLog{NULL};

    Vector3d xyz,lla,xyzR,llaR;//R:地球自传
    Vector3d vxyz,vxyzR;
    double tsdt,tsDrift;
    GnssTime ts0,ts;

    double elevationAngle,azimuthAngle;
    //for RTK:
    //vector<MSM4data*>里面是不同时刻的数据，0 always是最近的记录
    deque<MSM4data*> rtkData;
    deque<Measure*> measureDat;
    double prInterp[32], cpInterp[32];
    bool KalmanFirst{1};
    //电离层估计
    double Ii{0},IiP;
    char tip[128];
public:
    SV();
    SV(int id);
    ~SV();
    bool CheckEphemStates(GnssTime time,int ephemType);
    bool MeasureGood();
    bool ElevGood();
    bool IsMaskOn();
    bool RotateECEF(double tt);
    bool AddMmeasure(Measure *mesr);

    int CalcuelEvationAzimuth(Vector3d pos, Vector3d poslla);
    int CalcuTroposhphere(double elev,double azim);
    int CalcuInoshphere(double elev,double azim,Vector3d LLA,double time);
    int CorrectIT(Vector3d xyz,Vector3d lla,double time);
    void PrintInfo(int printType);
    void FPrintInfo(int printType);
    double InterpRtkData(double time, int sigInd);
    double InterpMeasere(int len,int power,int begin=0);
    double SmoothPr(int len,int begin=0);
    double SmoothKalman(int len,int begin=0);
    double SmoothKalman0();
    double DetectCycleSlip();

//private:
    double InterpLine(double xi,vector<double>&x,vector<double>&y);
};

class BeiDouSV:public SV{
//public:
//    //signal number:2,3,4;  8,9,10;  14,15,16;
//    SignalData B1_2I,B1_2Q,B1_2X,B3_6I,B3_6Q,B3_6X,B2_7I,B2_7Q,B2_7X;
public:
    BeiDouSV(int id);
    ~BeiDouSV();
    int DecodeSubFrame(uint32_t* dwrds){
        if(isBeiDouGEO)DecodeD2Frame1(dwrds);
        else DecodeD1(dwrds);
        return 1;
    }
    int DecodeD1(uint32_t* dwrds);
    int DecodeD2Frame1(uint32_t *dwrds);
//    SignalData* SignalTable(int index);

};

class GpsSV:public SV{
//public:
//    SignalData L1_1C,L1_1P,L1_1W,L2_2C,L2_2P,L2_2W,L2_2S,L2_2L,L2_2X,L5_5I,L5_5Q,L5_5X;
public:
    GpsSV(int id);
    ~GpsSV();
    int DecodeSubFrame(uint32_t* dwrds);
//    SignalData* SignalTable(int index);
};

class SvSys{
public:
    SysType type;
    vector<SV*> table;
    double tu;
    MatrixXd dataDebug;
    int xhead{0};
    char title[256];
    SvSys(SysType _type):type(_type){};
    void OpenClose(bool state){
        for(SV* sv:table)sv->open=state;
    }
    void sortElev(string tag){
        if("elev"==tag)
        sort(table.begin(),table.end(),[](SV* left,SV* right)->bool{ return left->elevationAngle>right->elevationAngle;});
    }
    void MakeDebug(int s){
        int Ni = table.size();
        dataDebug = MatrixXd(Ni,16);
        dataDebug.fill(0);
        xhead = 0;
        switch (s){
            case 3:
                for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->measureDat.front()->prMes;}
                xhead++;
                for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->measureDat.front()->cpMes;}
                xhead++;
            case 2:
               for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->I;}
                xhead++;
                for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->T;}
                xhead++;
            case 1:
                for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->elevationAngle*R2D;}
                xhead++;
                for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->azimuthAngle*R2D;}
                xhead++;
                for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->measureDat.front()->cno;}
                xhead++;
                for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->measureDat.front()->stdevPr;}
                xhead++;
                for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->measureDat.front()->stdevCp;}
                xhead++;
                for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->measureDat.front()->trkStat;}
                xhead++;
            case 0:
                for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->type;}
                xhead++;
                for(int i=0;i<Ni;i++){dataDebug(i,xhead)=table[i]->svId;}
                xhead++;
        }
    }
    void AddAnaData(VectorXd &v,int k=0){
        dataDebug.block(k,xhead++,v.rows(),1)=v;
    }
    void Show(){
        cout<<dataDebug<<endl;
    }

    template <typename T>
    VectorXd DiffDouble(T FUN,int sigId){
        int N = table.size();
        double fun0 = FUN(table[0],sigId);
        VectorXd result(N-1);
        for (int i = 1; i < N; ++i) {
            result(i-1) = fun0-FUN(table[i],sigId);
        }
        return result;
    }
    template <typename T>
    VectorXd DiffZero(T FUN,int sigId){
        int N = table.size();
        VectorXd result(N);
        for (int i = 0; i < N; ++i) {
            result(i) = FUN(table[i],sigId);
        }
        return result;
    }
    template <typename T>
    int SetValue(T FUN,int sigId,VectorXd &val){
        int N = table.size();
        for (int i = 0; i < N; ++i) {
            FUN(table[i],sigId) = val(i);
        }
        return 0;
    }
    MatrixXd GetE(Vector3d &pos,VectorXd &rs){
        int N = table.size();
        double r,t;
        MatrixXd result_T(3,N);
        Vector3d ei;
        for (int i = 0; i < N; ++i) {
            SV* sv = table[i];
            ei = sv->xyzR-pos;
            rs(i) = ei.norm();
            result_T.block<3,1>(0,i)<<ei/rs(i);
        }
        return result_T.transpose();
    }
    MatrixXd GetD(){
        int N = table.size();
        MatrixXd D(N-1,N);
        D.fill(0);
        for (int i =1; i <N; ++i) {
            D(i-1,i) =-1;
            D(i-1,0) = 1;
        }
        return D;
    }
};

class SvAll{
public:
    GNSS* gnss;
        Ionosphere ionKlob;

    vector<SvSys*> sysAll,sysUsed;
    vector<SV*> svAll, svUsed;

public:
    SvAll();
    void SetOpen(bool bds,bool gps);
//    SvAll(GNSS* gnss);
    ~SvAll();
    void InitAlloc();
    void UpdateEphemeris(char * subFrame);
    int UpdateSVs(string tag="NULL"){
       int count=0;
       vector<SV*> result;
       for(SvSys* sys:sysUsed){
           if("elev"==tag)sys->sortElev(tag);
           for(SV* sv:sys->table){
               result.push_back(sv);
               count++;
           }
       }
        svUsed.swap(result);
        return count;
    }
    //id begin from 1
    SV* GetSv(SysType type, int id){
        if(SYS_NULL==type)
            return nullptr;
        SvSys *sys =GetSys(type);
        if(id>sys->table.size()){
            printf("sv not found!\n");
            return nullptr;
        }
        return sys->table[id-1];
    }

    SvSys* GetSys(SysType type){
        SvSys* result= nullptr;
        for(SvSys* sys:sysAll){
            if(sys->type==type){
                result = sys;
                break;
            }
        }
        if(result== nullptr)
            printf("svSys Data not found!\n");
        return result;
    }
    SvSys* GetUsedSys(SysType type){
        SvSys* result = nullptr;
        for(SvSys* sys:sysUsed){
            if(sys->type==type)result = sys;
        }
        if(nullptr==result){
            sysUsed.push_back(new SvSys(type));
            result = sysUsed.back();
        }
        return result;
    }
    int AddToUsed(SV *sv);
private:
};


#endif