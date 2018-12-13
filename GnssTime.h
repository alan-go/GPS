//
// Created by root on 8/31/18.
//

#ifndef GPS_GNSSTIME_H
#define GPS_GNSSTIME_H

#include <string>
#include <time.h>
#include <cmath>
#include <string.h>
//#include "CommonInclude.h"
#define MAXLEAPS 64
const static double gpst0[]={1980,1, 6,0,0,0}; /* gps time reference */
const static double gst0 []={1999,8,22,0,0,0}; /* galileo system time reference */
const static double bdt0 []={2006,1, 1,0,0,0}; /* beidou time reference */
static double leaps[MAXLEAPS+1][7]={ /* leap seconds (y,m,d,h,m,s,utc-gpst) */
        {2017,1,1,0,0,0,-18},
        {2015,7,1,0,0,0,-17},
        {2012,7,1,0,0,0,-16},
        {2009,1,1,0,0,0,-15},
        {2006,1,1,0,0,0,-14},
        {1999,1,1,0,0,0,-13},
        {1997,7,1,0,0,0,-12},
        {1996,1,1,0,0,0,-11},
        {1994,7,1,0,0,0,-10},
        {1993,7,1,0,0,0, -9},
        {1992,7,1,0,0,0, -8},
        {1991,1,1,0,0,0, -7},
        {1990,1,1,0,0,0, -6},
        {1988,1,1,0,0,0, -5},
        {1985,7,1,0,0,0, -4},
        {1983,7,1,0,0,0, -3},
        {1982,7,1,0,0,0, -2},
        {1981,7,1,0,0,0, -1},
        {0}
};

class GnssTime {
public:
    time_t time{-1};
    double sec{0},tow{0},tod{0};
    int week{0};
//    static double leaps[MAXLEAPS+1][7]; /* leap seconds (y,m,d,h,m,s,utc-gpst) */

public:
    GnssTime();
    GnssTime(double tod);
    GnssTime(const double *ep);
    GnssTime(char *buffer, int len, bool utc2gps = 0);
    GnssTime(int week, double tow);
    void epoch2time(const double *ep);
    void epoch2time(const struct tm* ep);
    void time2epoch(double *ep);
    void utc2gpst();
    void time2tow();
    void gpst2utc();
    void Show(int mode=0);
    double operator-(const GnssTime &right);
    GnssTime operator+(double dt);
    GnssTime operator-(double dt);
    double operator+=(const double right);
    bool operator<(const GnssTime right);
    bool operator>(const GnssTime right);
};

#endif //GPS_GNSSTIME_H
