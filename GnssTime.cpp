//
// Created by root on 8/31/18.
//

#include "GnssTime.h"


GnssTime::GnssTime():time(-1),sec(-1) {}

GnssTime::GnssTime(const double *ep) {
    epoch2time(ep);
}

GnssTime::GnssTime(double tod):tod(tod+18.0){
}
GnssTime::GnssTime(int week, double tow):week(week),tow(tow) {
    epoch2time(gpst0);
    if (tow<-1E9||1E9<tow) tow=0.0;
    time+=86400*7*week+(int)tow;
    sec=tow-(int)tow;
}

GnssTime::GnssTime(char *head, int len, bool utc2gps):GnssTime() {
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
    if(utc2gps)utc2gpst();
}

void GnssTime::epoch2time(const struct tm*ep) {
    const int doy[]={1,32,60,91,121,152,182,213,244,274,305,335};
    int days,sec0,year=ep->tm_year+1900,mon=ep->tm_mon+1,day=ep->tm_mday;

    if (year<1970||2099<year||mon<1||12<mon) return;

    /* leap year if year%4==0 in 1901-2099 */
    days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3?1:0);
    sec0=ep->tm_sec;
    time=(time_t)days*86400+ep->tm_hour*3600+ep->tm_min*60+sec0;
    sec=0;
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

void GnssTime::time2tow(){
    GnssTime time0(gpst0);
    double secDiff = *this-time0;
    week = floor(secDiff/604800);
    tow = secDiff-week*604800;
    if (tow<-1E9||1E9<tow) tow=0.0;
}

void GnssTime::utc2gpst() {
    for (int i=0;leaps[i][0]>0;i++) {
        if ((*this-GnssTime(leaps[i]))>=0.0) {
            *this+=(-leaps[i][6]);
            return;
        }
    }
    time2tow();
}
void GnssTime::gpst2utc() {
    for (int i=0;leaps[i][0]>0;i++) {
        if ((*this-GnssTime(leaps[i]))>=0.0) {
            *this+=(leaps[i][6]);
            return;
        }
    }
}

double GnssTime::operator-(const GnssTime &right) {
    return difftime(time,right.time)+sec-right.sec;
}

GnssTime GnssTime::operator+(double dt) {
    GnssTime result = *this;
    result+=dt;
    return result;
}
GnssTime GnssTime::operator-(double dt) {
    GnssTime result = *this;
    result+=-dt;
    return result;
}
double GnssTime::operator+=(const double secAdd) {
    double tt;
    sec+=secAdd; tt=floor(sec); time+=(int)tt; sec-=tt;
    time2tow();
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

void GnssTime::time2epoch(double *ep) {
    const int mday[]={ /* # of days in a month */
            31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
            31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
    };
    int days,sec,mon,day;

    /* leap year if year%4==0 in 1901-2099 */
    days=(int)(time/86400);
    sec=(int)(time-(time_t)days*86400);
    for (day=days%1461,mon=0;mon<48;mon++) {
        if (day>=mday[mon]) day-=mday[mon]; else break;
    }
    ep[0]=1970+days/1461*4+mon/12; ep[1]=mon%12+1; ep[2]=day+1;
    ep[3]=sec/3600; ep[4]=sec%3600/60; ep[5]=sec%60+sec;
}
