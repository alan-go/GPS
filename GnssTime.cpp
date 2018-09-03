//
// Created by root on 8/31/18.
//

#include "GnssTime.h"


GnssTime::GnssTime():time(-1),sec(-1) {}

GnssTime::GnssTime(const double *ep) {
    epoch2time(ep);
}

GnssTime::GnssTime(int week, double tow) {
    epoch2time(gpst0);
    if (tow<-1E9||1E9<tow) tow=0.0;
    time+=86400*7*week+(int)tow;
    sec=tow-(int)tow;
}

GnssTime::GnssTime(char *head, int len):GnssTime() {
    double ep[6];
    char str[256],*p=str;

    if (strlen(head)<len-1) {
        return;
    }
    for (;*head&&--len>=0;) *p++=*head++; *p='\0';
    if (sscanf(str,"%lf %lf %lf %lf %lf %lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<6)
    {
        return;
    }
    if (ep[0]<100.0) ep[0]+=ep[0]<80.0?2000.0:1900.0;
    epoch2time(ep);
}

void GnssTime::epoch2time(const double *ep) {
    const int doy[]={1,32,60,91,121,152,182,213,244,274,305,335};
    int days,sec0,year=(int)ep[0],mon=(int)ep[1],day=(int)ep[2];

    if (year<1970||2099<year||mon<1||12<mon) return;

    /* leap year if year%4==0 in 1901-2099 */
    days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3?1:0);
    sec0=(int)floor(ep[5]);
    time=(time_t)days*86400+(int)ep[3]*3600+(int)ep[4]*60+sec0;
    sec=ep[5]-sec0;
}



void GnssTime::utc2gpst() {
    for (int i=0;leaps[i][0]>0;i++) {
        if ((*this-GnssTime(leaps[i]))>=0.0) {
            *this+=(-leaps[i][6]);
        }
    }
}
double GnssTime::operator-(const GnssTime &right) {
    return difftime(time,right.time)+sec-right.sec;
}

double GnssTime::operator+=(const double secAdd) {
    double tt;
    sec+=secAdd; tt=floor(sec); time+=(int)tt; sec-=tt;
    return double(time)+sec;
}
bool GnssTime::operator<(const GnssTime right) {
    if(*this-right<0)return true;
    return false;
}
bool GnssTime::operator>(const GnssTime right) {
    if(*this-right>0)return true;
    return false;
}
