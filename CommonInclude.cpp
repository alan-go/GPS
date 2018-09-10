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


int XYZ2LLA(Eigen::Vector3d XYZ, Eigen::Vector3d &LLA) {
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

void EarthRotate(Eigen::Vector3d in, Eigen::Vector3d &out, double dt){
    double omega = dt*Omega_e;
    Eigen::Matrix3d earthRotate(3,3);
    earthRotate<<cos(omega),sin(omega),0,-sin(omega),cos(omega),0,0,0,1;
    out = earthRotate*in;
}


