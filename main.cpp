#include "dm100.h"
#include "UBLOXM8L.h"
using namespace std;

int main()
{
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