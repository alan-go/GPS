//
// Created by alan on 18-12-13.
//
#include "GNSS.h"
int main(int argc, char* argv[]) {
    std::cout << "Hello, GNSS!" << std::endl;
    for(int i = 0;i<argc;i++)
        cout<<argv[i]<<endl;

    GNSS *gnss = new GNSS();
    //添加USB0(Ublox)
    gnss->AddSerial(0,0,"/dev/ttyUSB0",115200,1);
    gnss->GetSerial(0)->sendGGA=1;
    gnss->GetSerial(0)->ntripIn=1;

    //添加USB1()
    gnss->AddSerial(1,0,"/dev/ttyUSB1",115200,1);

    //log RTK数据
    gnss->rtkManager.logOpen=1;

    //Init参数:ephemType,usQianxun,Bds,Gps
    gnss->Init(0,1,1,1);

    //start
    gnss->StartGNSS();

    while('x'==getchar()) {
        cout << "stop capture." << endl;
        gnss->StopGNSS();
        sleep(1);
        break;
    }
    return 0;
}
