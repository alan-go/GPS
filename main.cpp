#include "dm100.h"
#include "UBLOXM8L.h"
using namespace std;

int main()
{
	cout<<atof("12.46l4,dji")<<endl;
	cout<<"start."<<endl;
	u_char a[8]={0x46,0xb6,0xf3,0xfd,0x87,0xb5,0x0b,0x41};
	/*printf("16: %02x,%02x,%02x,%02x,\n",*a,*(a+1),*(a+2),*(a+3));
	uint32_t aui = *(uint32_t *)a;
	printf("10: %d\n",aui);
	aui=aui>>19;
	aui=aui&0x7ff;
	printf("10: %d\n",aui);*/
	double tow  = *(double *)a;
	printf("%f",tow);


    UBLOXM8L ublox;
	//ublox.StartGPS("/dev/ttyUSB0");
	if('x'==getchar()){
	    cout<<"stop capture."<<endl;
	    stopUblox = 1;
	}
    //ublox.StopGPS();

    return 0;
}