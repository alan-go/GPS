#include "dm100.h"
#include "UBLOXM8L.h"
using namespace std;

inline uint32_t Read1Word(uint32_t word, int head, int length) {
	uint32_t result = (word<<head)>>(32-length);
	return result;
}


inline uint32_t Read2Word(uint32_t* word, int head0,int head1, int length0,int length1) {
	uint32_t bit32x2[2];
	uint64_t result;
	bit32x2[0] = word[0]<<head0>>(32-length0)<<length1;
	bit32x2[1] = word[1]<<head1>>(32-length1);
	result = bit32x2[0]|bit32x2[1];
	printf("read2:%08x %08x -> %08x\n",bit32x2[0],bit32x2[1],result);
	return result;
}


int main()
{
	cout<<atof("12.46l4,dji")<<endl;
	cout<<"start."<<endl;
	//u_char a[80]={0x79,0x13,0x90,0x38,0x87,0x40,0xf8,0x1b, 0xD9 ,0xE6,0x50 ,0x00 ,0xBF ,0x23 ,0xF8 ,0x3A};
	u_char a[80]={0x79,0x13,0x90,0x38,0x87,0x40,0xf8,0x1b, 0xD9 ,0xE6,0x50 ,0x00 ,0xBF ,0x23 ,0xF8 ,0x3A};
	//1 0111 0011 1111 0
	/*printf("16: %02x,%02x,%02x,%02x,\n",*a,*(a+1),*(a+2),*(a+3));
	uint32_t aui = *(uint32_t *)a;
	printf("10: %d\n",aui);
	aui=aui>>19;
	aui=aui&0x7ff;
	printf("10: %d\n",aui);*/
	double tow  = *(double *)a;
	printf("%f",tow);

	uint32_t bit32 = *(uint32_t*)a;
	uint32_t frameId = Read1Word(bit32,17,3);
	printf("\nbit32:%08x,Frame:%u",bit32,frameId);

	uint32_t *bit32x2 = (uint32_t *)a;
	uint32_t word[3],head[3],length[3];
	printf("\nbit32:%08x,%08x\n",bit32x2[2],bit32x2[3]);
	uint64_t bit64 = Read2Word(bit32x2+2,15,2,9,8);
	printf("bit64 = %llu\n",bit64);

	uint32_t ii32 = -23;
	int32_t ii = -23;
	double dd32 = ((int32_t)ii32)*pow(2,-24);
	printf("dd32 =  %u,%d, %lf,%lf",ii32, ii32, dd32,(double)ii);



    UBLOXM8L ublox;
	//ublox.StartGPS("/dev/ttyUSB0");
	if('x'==getchar()){
	    cout<<"stop capture."<<endl;
	    stopUblox = 1;
	}

    return 0;
}