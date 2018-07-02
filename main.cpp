#include "dm100.h"
#include "UBLOXM8L.h"
using namespace std;

int main()
{
    Eigen::MatrixXd earthRotate(3,3);
    earthRotate<<1,2,3,2,2,3,2,7,6;
    Eigen::Vector3d xyz3;
    Eigen::Vector4d xyz4;
    xyz3<<1,2,3;
    xyz4<<2,1,3,3;
    Eigen::Vector3d xyz43 = xyz4.head(3);
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
    SvInfo *aa;
    aa = new SvInfo;
    cout<<"sqa="<<sq_M_miu<<"  Omegae"<<Omega_e<<endl;
    delete(aa);

    UBLOXM8L ublox;
	ublox.StartGPS("/dev/ttyUSB0");
	if('x'==getchar()){
	    cout<<"stop capture."<<endl;
	    stopUblox = 1;
	}

    return 0;
}