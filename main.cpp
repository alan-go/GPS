#include "dm100.h"
#include "UBLOXM8L.h"
using namespace std;

int main()
{
	cout<<atof("12.46l4,dji")<<endl;
	cout<<"start."<<endl;

    UBLOXM8L ublox;
	//ublox.StartGPS("/dev/ttyUSB0");
	if('x'==getchar()){
	    cout<<"stop capture."<<endl;
	    stopUblox = 1;
	}

    return 0;
}