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
    for(i=0;i<=n-1;i++)
    { a[i]=y[i];
        for(j=0;j<=n-1;j++)
            if(j!=i) a[i]*=(xx-x[j])/(x[i]-x[j]);
        yy+=a[i];
    }
    free(a);
    return yy;
}


