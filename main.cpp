#include "dm100.h"
#include "UBLOXM8L.h"
using namespace std;

int main()
{
	char t[8];
	t[0] = 0xb5;
	t[1] = 'm';
	printf("%x\t%x\n",(unsigned char)t[0],t[1]);
	if((u_char)t[0]==0xb5){
		cout<<"((unsigned char)t[0]==0xb5)"<<endl;
	}
	cout<<atof("12.46l4,dji")<<endl;
	cout<<"start"<<endl;
	UBLOXM8L ublox;
	ublox.StartGPS("/dev/ttyUSB0");
	if('x'==getchar()){
	    cout<<"stop capture."<<endl;
	    stopUblox = 1;
	}
    //ublox.StopGPS();

    return 0;
}