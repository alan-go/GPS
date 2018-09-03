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