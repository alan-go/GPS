//
// Created by root on 8/1/18.
//

#include "GNSS.h"
#include <stdio.h>
using namespace std;
using namespace Eigen;



int main(){
    BeiDouSV sv;
    vector<double >xx,yy;
    xx.push_back(553185);
    xx.push_back(553184);
    xx.push_back(553183);
    xx.push_back(553182);
    xx.push_back(553181);

    yy.push_back(38231393.730568);
    yy.push_back(38231389.495611);
    yy.push_back(38231385.242785);
    yy.push_back(38231380.989959);
    yy.push_back(38231376.737132);

    double t=553185.002000,pr=38230924.070819;
    double yt = sv.InterpLine(t,xx,yy);
    printf("yt = %f\n",yt);


    printf("\nin testing\n");
    GNSS gnss;
//    gnss.useGPS = false;
    gnss.useBeiDou = false;
//    gnss.useQianXun = false;

    u_char bbb[36] = {0xD3 ,0x00 ,0x13 ,0x3E ,0xD7 ,0xD3 ,0x02 ,0x02 ,0x98 ,0x0E ,0xDE ,0xEF ,0x34 ,0xB4 ,0xBD ,0x62
            ,0xAC ,0x09 ,0x41 ,0x98 ,0x6F ,0x33 ,0x36 ,0x0B ,0x98 ,0xd3, 0x00};
    char aaa[36];
    memcpy(aaa,bbb,30);

    u_char ttt[256] =     { 0xD3 ,0x00 ,0x8A ,0x43 ,0x20 ,0x33 ,0x5B ,0x36 ,0x96 ,0x02 ,0x00 ,0x00 ,0x00 ,0x0B ,0x26 ,
                            0xC0 ,0x00 ,0x00 ,0x00 ,0x00 ,0x20 ,0x20 ,0x00 ,0x00 ,0x7F ,0xFF ,0xA4 ,0xA6 ,0xA5 ,0xA1 ,
                            0xA5 ,0x24 ,0xA3 ,0x26 ,0x9F ,0x03 ,0xDF ,0xB3 ,0x2D ,0x77 ,0xE1 ,0x90 ,0x52 ,0x78 ,0x99 ,
                            0x56 ,0x36 ,0x0E ,0x1E ,0xCC ,0x4A ,0x66 ,0x59 ,0xBC ,0xE7 ,0x9A ,0xAC ,0xB6 ,0x0C ,0x1E ,
                            0xE5 ,0x41 ,0x82 ,0x61 ,0xA0 ,0xCB ,0xC7 ,0x5A ,0xBE ,0xF5 ,0xFB ,0x57 ,0xB8 ,0x7E ,0x89 ,
                            0x51 ,0x90 ,0x2C ,0xA8 ,0x30 ,0x48 ,0xF3 ,0x01 ,0x3D ,0xD9 ,0x7D ,0x44 ,0x44 ,0x03 ,0xC1 ,
                            0x2F ,0xC1 ,0x5F ,0xDF ,0x2D ,0xF7 ,0x08 ,0x45 ,0xD2 ,0x24 ,0x27 ,0x80 ,0x73 ,0x4E ,0xC1 ,
                            0xFC ,0x7E ,0x01 ,0x4C ,0xC6 ,0x15 ,0x6A ,0x57 ,0xD8 ,0x28 ,0xA0 ,0x47 ,0x32 ,0xFF ,0xFF ,
                            0xFF ,0xFF ,0xFF ,0xFF ,0xFF ,0xFF ,0x80 ,0x00 ,0x63 ,0x05 ,0x6C ,0x5E ,0xFE ,0x74 ,0x58 ,
                            0xCE ,0x50 ,0xE7 ,0x35 ,0x8B ,0x00 ,0xc4 ,0xce ,0xf3};

    u_char ttt2[256] = {0xd3, 0x00, 0xa7, 0x43, 0x28, 0xf4, 0x40, 0x28, 0x51, 0x00, 0x00, 0x40, 0xd1, 0x90, 0x62, 0x0a,
                        0x00, 0x00, 0x00, 0x00, 0x20, 0x40, 0x00, 0x00, 0x7f, 0xff, 0xfa, 0x1a, 0x72, 0x62, 0x72, 0x22,
                        0x62, 0x32, 0x62, 0x2a, 0x43, 0x9e, 0x4a, 0x60, 0xb8, 0xd6, 0xb8, 0x36, 0xd0, 0xc6, 0x9d, 0xed,
                        0xb3, 0xe2, 0x1a, 0xc2, 0xfe, 0x41, 0x5c, 0x7f, 0x30, 0xf4, 0x11, 0xc2, 0x66, 0x42, 0x8c, 0x43,
                        0x66, 0x4b, 0xc9, 0xb1, 0x82, 0x22, 0xfa, 0xb6, 0xb6, 0x4d, 0x48, 0xbc, 0x8c, 0xf8, 0x47, 0x0b,
                        0xb3, 0x14, 0x68, 0x3f, 0xf0, 0x7a, 0xd2, 0x71, 0x6d, 0x89, 0xa6, 0xa7, 0x48, 0x5f, 0x33, 0x24,
                        0x33, 0xb2, 0xff, 0x92, 0x0b, 0xe5, 0x37, 0x24, 0x3c, 0x0a, 0x90, 0x69, 0x6c, 0x41, 0xfb, 0xb1,
                        0x0f, 0x00, 0x6e, 0xd5, 0xc3, 0xbb, 0x67, 0x7e, 0x49, 0xfa, 0xd9, 0x2d, 0x94, 0x6c, 0x81, 0xd1,
                        0xb5, 0x16, 0x76, 0x15, 0xae, 0xf8, 0x5c, 0x71, 0x8d, 0x17, 0x6a, 0x33, 0xf8, 0xbf, 0xff, 0xff,
                        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xf8, 0x00, 0x00, 0x67, 0x95, 0x34, 0xdd, 0x5d, 0x74,
                        0xe5, 0x7e, 0x37, 0x61, 0x75, 0xb5, 0xe3, 0x7e, 0x57, 0x80, 0x66, 0xf4, 0x6c };
    char msm1074[256];
    memcpy(msm1074,ttt,256);


//    gnss.rtkManager.TestParase(msm1074,256);
//
//    gnss.rtkManager.TestParase(aaa,30);
//    cout<<"ECEF = "<<gnss.rtkManager.ECEF_XYZ<<endl;
//    printf("Y = %lf\n",gnss.rtkManager.ECEF_XYZ(1));

    string serialPort = "/dev/ttyUSB0";
    gnss.StartGNSS(serialPort,115200);
    getchar();

    return 0;
}
