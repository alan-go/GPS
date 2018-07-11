#include "dm100.h"
#include "UBLOXM8L.h"
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
    SvInfo sv;
    sv.orbit.toe = 244800;
    sv.orbit.sq_a = 5153.65531;
    sv.orbit.e = 0.005912038265;
    sv.orbit.i0 = 0.9848407943;
    sv.orbit.Omega0 = 1.038062244;
    sv.orbit.omega = -1.717457876;
    sv.orbit.M0 = -1.064739758;
    sv.orbit.dtn = 4.249105564e-9;
    sv.orbit.IDOT = 7.422851197e-51;
    sv.orbit.OmegaDot = -8.151768125e-9;
    sv.orbit.Cuc = 3.054738045e-7;
    sv.orbit.Cus = 2.237036824e-6;
    sv.orbit.Crc = 350.53125;
    sv.orbit.Crs = 2.53125;
    sv.orbit.Cic = -8.381903172e-8;
    sv.orbit.Cis = 8.940696716e-8;
    sv.CalcuECEF(239050.7223);
    sv.PrintInfo(1);
    Vector3d xyz,LLA;
    xyz<<-2267521,5008960,3221750;
    UbloxSolver sssver;
    sssver.XYZ2LLA(xyz,LLA);
    cout<<"xyz"<<xyz<<endl;
    cout<<"LLA="<<(180/GPS_PI)*LLA<<endl;

    double a = 34;
    Vector2d v2d(a,54);
    cout<<v2d<<endl;
    MatrixXd earthRotate(3,3);
    earthRotate<<1,2,3,2,2,3,2,7,6;
    Vector3d xyz3;
    Vector4d xyz4;
    xyz3<<1,2,3;
    xyz4<<2,1,3,3;
    Vector3d xyz43 = xyz4.head(3);
    double r = (xyz3 - xyz43).squaredNorm();
    cout<<xyz3<<xyz43<<r<<endl;

    cout <<earthRotate<<endl;
    SvInfo infoa[2],infob[3];
    infoa[0].a0 =1.98;
    infob[2] = infoa[0];
    infoa[0].a0 = 235;
    cout<<"infoa/b.a0 = "<<infoa[0].a0<<infob[2].a0<<endl;

	cout<<atof("12.46l4,dji")<<endl;
	cout<<"start."<<endl;






    UBLOXM8L ublox;
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
            ublox.ScanSerialData(dat,128);
        }
        inF.close();

    } else if('r'==command){
        ublox.StartGPS(serialPort);
        if('x'==getchar()){
            cout<<"stop capture."<<endl;
            stopUblox = 1;
        }
    }


    return 0;
}