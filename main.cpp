
#include "GNSS.h"
using namespace std;
using namespace Eigen;

//string serialPort = "/dev/tty.usbserial";
string serialPort = "/dev/ttyUSB0";
std::ofstream outF;
bool stopLog = 0;
char saveDataName[64];

int MakeGGA(char *gga, Vector3d lla, GnssTime gpsTime) {
    GnssTime ggaTime = gpsTime;
    double h=0,ep[6],dms1[3],dms2[3],dop=1.0;
    int solq = 1;
    char *p=gga,*q,sum;

    ggaTime.gpst2utc();
    if (ggaTime.sec>=0.995) {ggaTime.time++; ggaTime.sec=0.0;}
    ggaTime.time2epoch(ep);

    deg2dms(fabs(lla(0))*180.0/GPS_PI,dms1);
    deg2dms(fabs(lla(1))*180.0/GPS_PI,dms2);
    p+=sprintf(p,"$GPGGA,%02.0f%02.0f%05.2f,%02.0f%010.7f,%s,%03.0f%010.7f,%s,%d,%02d,%.1f,%.3f,M,%.3f,M,%.1f,",
               ep[3],ep[4],ep[5],dms1[0],dms1[1]+dms1[2]/60.0,lla(0)>=0?"N":"S",
               dms2[0],dms2[1]+dms2[2]/60.0,lla(1)>=0?"E":"W",solq,
               6,dop,lla(2),h,1);
    for (q=(char *)gga+1,sum=0;*q;q++) sum^=*q; /* check-sum */
    p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
    return p-(char *)gga;
}

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

void test(GNSS *gnss){
    GnssTime timeg(2018,208953.253);
    char gganow[128];
    MakeGGA(gganow,Eigen::Vector3d(0.69777,2.03050,52.91082),timeg);
    printf("%s\n",gganow);

    char a[12] = "1.36e1";
    printf("a= %f\n", (str2num(a, 6)));

    MatrixXf D(10,11),temp;

    VectorXf R_(11);
    D.fill(0);
    for (int i = 0; i < 10; ++i) {
        D(i,i+1) = -1;
        R_(i) = i+2;
    }
    D.block<10,1>(0,0).fill(1);
    MatrixXf R(R_.asDiagonal());
    R.colwise();
    cout<<R<<endl;
    cout<<R_<<endl;

    cout<<D<<endl;
    cout<<D*R*D.transpose()<<endl;
    MatrixXd B(4,3);
    Vector3d Bb{3,2.3,8};
    VectorXd B6(6);
    B6.head(6)<<Bb,Bb;
    cout<<"B6\n"<<B6<<endl;

    int foo[5]{1,2,3,4,5};
    int* pfoo = &foo[3];
    cout<<foo <<"pfoo"<<pfoo<<endl;
    *pfoo = 464;
    printf("f3 %d\n", foo[3]);
    MatrixXd A = MatrixXd::Random(6,3);
    VectorXd yp(6);yp<<2,5,3,2,6,7;
    cout <<"ypa"<<A.colPivHouseholderQr().solve(yp)<<endl;

    VectorXd xii(5);
    xii.fill(3);


}

int main()
{
    GNSS *gnss = new GNSS();
    gnss->AddSerial(0,0,"/dev/ttyUSB0",115200,true,false);
    gnss->GetSerial(0)->sendGGA=1;
    gnss->AddSerial(1,1,"/dev/ttyUSB1",115200,true,false);
    gnss->rtkManager.logOpen=1;
//    test(gnss);
//    return 0;
    for (int i = 0; i < 5; ++i) gnss->svMaskBds[i]=0;
    gnss->log = fopen("../log/log.txt","w");
    gnss->logDebug = fopen("../log/logDebug.txt","w");
    gnss->useQianXun = true;
    //ephem,qianxun,bds,gps
    gnss->Init(0,1,1,1);
//    gnss->Init(0,0,1,0);
//    gnss->Init(0,0,0,1);


    printf("command:\nw : write serial.\nd : from data.\nr : from receiver.\n");
    char command = getchar();
//    char command = 'd';
    if('d'==command) {
        FILE *fp;
        char name[128],dat[512],temp[256],tempc;
        string ss = "0921_13_02";
//        string ss = "0802-1";
//        string ss = "0708-2.data";
//        string ss = "0823";
//        string ss = "0815-2";//soho novatal
//        scanf("%s",name);

        string ssData = "../data/" + ss + ".data";
        string ssRTK = "../data/" + ss + ".rtk";
        printf("open file name:%s\n",ssRTK.data());
        if(fp = fopen(ssRTK.data(),"rb")){
        int k = 0;
        while (!feof(fp)){
            dat[k++] = fgetc(fp);
            if(0xd3==(u_char)dat[k-2]&&0==(dat[k-1]>>2)){
                gnss->rtkManager.ParaseRTK(dat,k);
                k=2;
                dat[0] = dat[k-2];
                dat[1] = dat[k-1];
                memset(dat+2,0,510);
            }
        }
        fclose(fp);
        } else printf("open rtk data failed \n");

        if(fp = fopen(ssData.data(),"rb")){
        while (128==fread(dat,1,128,fp)){
            gnss->GetSerial(0)->ScanSerialData(dat,128);
        }
        fclose(fp);
        } else printf("open data failed \n");

    } else if('r'==command){
        gnss->StartGNSS();
        char c = getchar();
        if('x'==getchar()){
            cout<<"stop capture."<<endl;
        }
    }else if('w'==command){
        boost::asio::serial_port* sp = gnss->GetSerial(0)->sp_ ;
        char bssUbx[64]="BSS-UBX-115200-5HZ-PVT-RHXZ-D1/r/n";
        boost::asio::io_service ios;
        sp = new boost::asio::serial_port(ios, serialPort);
        sp->set_option ( boost::asio::serial_port::baud_rate ( 115200 ) );
//        sp->write_some(bssUbx);
        boost::asio::write(*sp, boost::asio::buffer(bssUbx, strlen(bssUbx)+2));

    }
    fclose(gnss->log);
    sleep(2);
    cout<<"Quit?"<<endl;
    sleep(3);
//    getchar();
    return 0;
}