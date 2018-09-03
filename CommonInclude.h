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

//todo test this
namespace TOOL{

}

extern double str2num(char *head, int len);

extern double lagrange(double *x,double *y,double xx,int n);     /*拉格朗日插值算法*/

#endif //GPS_COMMONINCLUDE_H
