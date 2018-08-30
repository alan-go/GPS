
#include "GNSS.h"
using namespace std;
using namespace Eigen;

//string serialPort = "/dev/tty.usbserial";
string serialPort = "/dev/ttyUSB0";
std::ofstream outF;
bool stopLog = 0;
char saveDataName[64];

void *LogData(void *fileName){
    outF.open((char*)fileName,std::ofstream::binary);
    try {
        boost::asio::io_service ios;
        boost::asio::serial_port sp ( ios, serialPort );
        sp.set_option ( boost::asio::serial_port::baud_rate ( 115200 ) );
        printf ( "successfully opened port %s\n", serialPort.c_str() );

        while ( !stopLog ) {
            char tmp[256];
            auto transferred = sp.read_some ( boost::asio::buffer ( tmp ) );
            if ( transferred <= 0 ) {
                usleep ( 100 );
                printf ( "serial port return 0\n" );
                continue;
            }

            for(int k = 0;k<transferred;k++){
                printf("%02x ",(u_char) tmp[k]);
            }
            printf("\n");
//            printf("write%d\n",transferred);
            outF.write(tmp,transferred);
        }
        sp.close();
        outF.close();
    } catch ( ... ) {
        printf ( "failed to open serial port\n" );
    }
}

int main()
{
    GNSS *gnss = new GNSS();
    gnss->log = fopen("../log/log.txt","w");
//    gnss->log = fopen("../log/log0829-b1.txt","w");


    Vector3d pc[5],LLA[5];
    pc[0]<<-32353.678517,  27026.346058,  -1095.094873 ;
    pc[1]<<  4354.572550,  41921.050023,   1042.240678;
    pc[2]<< -14749.147484,  39509.822234,   -650.340449 ;
    pc[3]<< -39599.198043,  14448.060190,   -564.747852 ;
    pc[4]<<  21928.236952,  36015.925416,    895.460497 ;

//    PC01 -32353.678517,  27026.346058,  -1095.094873    312.878333
//    PC02   4354.572550,  41921.050023,   1042.240678    -16.403391
//    PC03 -14749.147484,  39509.822234,   -650.340449   -299.578701
//    PC04 -39599.198043,  14448.060190,   -564.747852   -273.299420
//    PC05  21928.236952,  36015.925416,    895.460497    315.537101
//    PC06  -6792.119072  30515.529479  28725.661181     89.249118
//    PC07 -22186.980236  22343.659926 -28133.021383     76.004768

    SVs svs; NtripRTK *rtk;
    PosSolver svlr(svs,rtk,gnss);
    for (int i=0;i<5;i++){
        svlr.XYZ2LLA(pc[i]*1000,LLA[i]);
        printf("\nBDS_GEO_LLA, %02d\n%lf\n%lf\n%lf\n",i+1,LLA[i](0)*180/GPS_PI,LLA[i](1)*180/GPS_PI,LLA[i](2));
    }



//    gnss->useGPS = false;
//    gnss->useBeiDou = false;
    gnss->useQianXun = false;

//    gnss->StartGNSS("null",115200);

    printf("command:\nl : log data.\nd : from data.\nr : from receiver.\n");
    char command = getchar();
//    char command = 'd';
    if('l'==command){
        printf("start write file %s",saveDataName);
        pthread_t logThread = 0;
        pthread_create(&logThread, nullptr,LogData,saveDataName);

        while (command=getchar()){
            if('x'==command){
                stopLog = 1;
                pthread_join(logThread, nullptr);
                printf("stopLogging");
                break;
            }
        }

    } else if('d'==command) {
        ifstream inF;
        char name[128],dat[128];
        printf("open file name:");
//        string ss = "../data/0802-1";
//        string ss = "../data/0708-2";
//        string ss = "../data/0823";
        string ss = "../data/0815-2";//soho novatal
//        scanf("%s",name);
        inF.open(ss, std::ifstream::binary);
        while (!inF.eof()){
            inF.read(dat,128);
//            sleep(1);
            gnss->serialDataManager.ScanSerialData(dat,128);
        }
        inF.close();

    } else if('r'==command){
        gnss->StartGNSS(serialPort,115200);
        if('x'==getchar()){
            cout<<"stop capture."<<endl;
//            gnss.stop;
//            stopUblox = 1;
        }
    }

    fclose(gnss->log);
    sleep(2);
    cout<<"Quit?"<<endl;
//    getchar();
    return 0;
}