//
// Created by alan on 18-12-13.
//
#include "GNSS.h"
int main() {
    std::cout << "Hello, GNSS!" << std::endl;
    double ep[6]={2018,11,13,9,9,0};
    GnssTime time(ep);
    time.utc2gpst();
    time.Show();
    return 0;
}
