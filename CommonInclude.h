//
// Created by root on 8/31/18.
//

#ifndef GPS_COMMONINCLUDE_H
#define GPS_COMMONINCLUDE_H

#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <mutex>
#include <queue>
#include <chrono>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#define NGPS 64
#define NBeiDou 64
//PI ???GPs_PI
constexpr static double GPS_PI = 3.1415926535898;
constexpr static double M_miu = 3.986004418e14;
constexpr static double Omega_e = 7.2921150e-005;//SYS_GPS?Beidou
constexpr static double Light_speed = 299792358.0;
constexpr static double Earth_a = 6378137.0;  //地球长半轴
constexpr static double Earth_f = 3.352810664747481e-003;   //基准椭球体的极扁率  f = 1/298.257223563
constexpr static double Earth_ee = 6.694379990141317e-003;   //偏心率e   e^2 = f(2-f)


extern double str2num(char *head, int len);

extern double lagrange(double *x,double *y,double xx,int n);     /*拉格朗日插值算法*/

extern int XYZ2LLA(Eigen::Vector3d XYZ, Eigen::Vector3d &LLA);

extern void deg2dms(double deg, double *dms);
#endif //GPS_COMMONINCLUDE_H
