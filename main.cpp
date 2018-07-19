
#include "GNSS.h"
using namespace std;
using namespace Eigen;

//string serialPort = "/dev/tty.usbserial";
string serialPort = "/dev/ttyUSB0";
std::ofstream outF;
bool stopLog = 0;

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
    GNSS gnss;

    printf("command:\nl : log data.\nd : from data.\nr : from receiver.\n");
//    char command = getchar();
    char command = 'd';
    if('l'==command){
        printf("\nfile name : ");
        char name[64];
        scanf("%s",name);
        getchar();
        printf("start write file %s",name);
        pthread_t logThread = 0;
        pthread_create(&logThread, nullptr,LogData,name);

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
        string ss = "0708-3";
//        scanf("%s",name);
        inF.open(ss, std::ifstream::binary);
        while (!inF.eof()){
            inF.read(dat,128);
            gnss.serialDataManager.ScanSerialData(dat,128);
        }
        inF.close();

    } else if('r'==command){
        gnss.StartGNSS(serialPort,115200);
        if('x'==getchar()){
            cout<<"stop capture."<<endl;
//            gnss.stop;
//            stopUblox = 1;
        }
    }


    return 0;
}