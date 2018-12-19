//
// Created by root on 9/1/18.
//
#include "CommonInclude.h"

double str2num(char *head, int len){
    double value;
    char str[256],*p=str;
    for (;*head&&--len>=0;head++) *p++=*head=='d'||*head=='D'?'E':*head; *p='\0';
    return sscanf(str,"%lf",&value)==1?value:0.0;
}

double lagrange(double *x,double *y,double xx,int n){
    int i,j;
    double *a,yy=0.0;    /*a作为临时变量，记录拉格朗日插值多项式*/
    a=(double *)malloc(n*sizeof(double));
    for(i=0;i<=n-1;i++) {
        a[i]=y[i];
        for(j=0;j<=n-1;j++)
            if(j!=i) a[i]*=(xx-x[j])/(x[i]-x[j]);
        yy+=a[i];
    }
    free(a);
    return yy;
}
double lineIntp(double *x,double *y,double xx,int n){

}

int XYZ2LLA(Eigen::Vector3d &XYZ, Eigen::Vector3d &LLA) {
    double x,y,z,r,sinA,cosA,sinB,cosB,N,h,lati;
    int i;
    double v_xyz[3],v_enu[3],a_xyz[3],a_enu[3];
    double lati_f,long_f;   // N_geoi表中的经纬度，以角度记，均加到正值处以10
    double lati_w,long_w;   //权重
    int lati_index,long_index,long_index_n;  //用于N_geoid中的下标
    //input
    x = XYZ(0);
    y = XYZ(1);
    z = XYZ(2);
    v_xyz[0] = v_xyz[0];
    v_xyz[1] = v_xyz[1];
    v_xyz[2] = v_xyz[2];
    a_xyz[0] = a_xyz[0];
    a_xyz[1] = a_xyz[1];
    a_xyz[2] = a_xyz[2];
    //起始迭代点
    r = sqrt(x*x+y*y);
    h = 0;
    lati = 0;
    for (i=0;i<5;i++){
        sinA = sin(lati);
        cosA = cos(lati);
        N = Earth_a / sqrt(1-sinA*sinA*Earth_ee);
        h = r/cosA - N;
        lati = atan(z/(r*(1-Earth_ee*N/(N+h))));
//        printf("lati = %lf\n",lati);
    }
    //output Longtitude Latitude High  SS[8] 转换矩阵
    //0,1,2:纬,经,高
    LLA(0) = lati;
    LLA(1) = atan2(y,x);
    LLA(2) = h;

    //给转换矩阵赋值
    sinB = y/r;
    cosB = x/r;
    double SS[8];
    SS[0] = -sinB;
    SS[1] = cosB;
    SS[2] = 0;
    SS[3] = -sinA*cosB;
    SS[4] = -sinA*sinB;
    SS[5] = cosA;
    SS[6] = cosA*cosB;
    SS[7] = cosA*sinB;
    SS[8] = sinA;

    //XYZ向速度，加速度ENU 之后有时间自己实现吧，下面就是一个矩阵乘法然后赋值
    //SS*v_xyz = v_enu   SS*a_xyz = a_enu
    return 0;
}

void deg2dms(double deg, double *dms)
{
    double sign=deg<0.0?-1.0:1.0,a=fabs(deg);
    dms[0]=floor(a); a=(a-dms[0])*60.0;
    dms[1]=floor(a); a=(a-dms[1])*60.0;
    dms[2]=a; dms[0]*=sign;
}

/* convert deg-min-sec to degree -----------------------------------------------
* convert degree-minute-second to degree
* args   : double *dms      I   degree-minute-second {deg,min,sec}
* return : degree
*-----------------------------------------------------------------------------*/
double dms2deg(const double *dms)
{
    double sign=dms[0]<0.0?-1.0:1.0;
    return sign*(fabs(dms[0])+dms[1]/60.0+dms[2]/3600.0);
}

void EarthRotate(Eigen::Vector3d in, Eigen::Vector3d &out, double dt){
    double omega = dt*Omega_e;
    Eigen::Matrix3d earthRotate(3,3);
    earthRotate<<cos(omega),sin(omega),0,-sin(omega),cos(omega),0,0,0,1;
    out = earthRotate*in;
}
void ShowV3(Eigen::Vector3d v3,char *tip){
    printf("^^ %s,,,%.10f,%.10f,%.10f,%.10f\n", tip,v3(0),v3(1),v3(2),v3.norm());
}

double GetFreq(SysType type, int sigInd, bool lambda){
    int i = sigInd+1;
    double result = -1;
    switch(type){
        case SYS_BDS:
            if(i>=2&&i<=4)result = 1561.098;//B1
            if(i>=8&&i<=10)result = 1268.52;//B3
            if(i>=14&&i<=16)result = 1207.14;//B2
            break;
        case SYS_GPS:
            if((i>=2&&i<=4)||(i>=30&&i<=32))result = 1575.42;//L1
            if((i>=8&&i<=10)||(i>=15&&i<=17))result = 1227.6;//L2
            if(i>=22&&i<=24)result = 1176.45;//L5
            break;
        default:
            return -1;
    }
    result*=1e6;
    result = lambda?(Light_speed/result):result;
    return result;
}
void LLA2XYZ(const Eigen::Vector3d &lla,Eigen::Vector3d &xyz )
{
    double sinp=sin(lla(0)),cosp=cos(lla(0)),sinl=sin(lla(1)),cosl=cos(lla(1));
    double e2=FE_WGS84*(2.0-FE_WGS84),v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);

    xyz(0)=(v+lla(2))*cosp*cosl;
    xyz(1)=(v+lla(2))*cosp*sinl;
    xyz(2)=(v*(1.0-e2)+lla(2))*sinp;
}

void XYZ2ENU(const Eigen::Vector3d &xyz,const Eigen::Vector3d &lla,Eigen::Vector3d &enu){
    double sinp=sin(lla[0]),cosp=cos(lla[0]),sinl=sin(lla[1]),cosl=cos(lla[1]);
    Eigen::Matrix3d SS;
    SS<<-sinl,      cosl,   0,
    -sinp*cosl, -sinp*sinl, cosp,
    cosp*cosl,  cosp*sinl,  sinp;
    enu = SS*xyz;
}
Solution FindSol(SolutionDeque &sols, double t,double dt,std::string tag )
{
    if(tag=="tod"){
        for(Solution sol:sols){
            if(abs(sol.time.tod-t)<dt)
                return sol;
        }
    }
}

//
//template <typename T>
//int FindTimeDate(T list,GnssTime &time){
//    printf("tempFindData \n");
//    /* binary search ->find head*/
//    int N = list.size();
//    int i=0,k=0,j;
//    for (i=0,j=N-1;i<j;) {
//        k=(i+j)/2;
//        if (list[k].time<time) i=k+1; else j=k;
//    }
//    if(i<1||i>N)
//        return -1;
//    return i;
//}
