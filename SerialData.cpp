//
// Created by root on 7/13/18.
//

#include "SerialData.h"
#include "GNSS.h"

SerialData::SerialData() :
flag(0),count(0),lengthNMEA(0),lengthUBX(0),baudRate(115200),stopCapture(false),showData(false){
    sp_ = nullptr;
}

SerialData::~SerialData() {

}

void SerialData::StartCapture(const std::string serialPort, unsigned int baudRate) {
    try {
        boost::asio::io_service ios;
        sp_ = new boost::asio::serial_port(ios, serialPort);

        sp_->set_option ( boost::asio::serial_port::baud_rate ( baudRate ) );
        printf ( "successfully opened port %s\n", serialPort.c_str() );

        while ( !stopCapture ) {
            char tmp[256];
            auto transferred = sp_->read_some ( boost::asio::buffer ( tmp ) );
            if ( transferred <= 0 ) {
                usleep ( 100 );
                printf ( "serial port return 0\n" );
                continue;
            }
//                printf ( "transferred = %d , flag = %d\n",transferred, flag );
            ScanSerialData(tmp,transferred);
        }
        sp_->close();
    } catch ( ... ) {
        printf ( "failed to open serial port\n" );
        return;
    }
}

void SerialData::ScanSerialData(char *tmp, int transferred) {
    for ( int i = 0; i<transferred; i++ ) {
        if ( '$'==tmp[i]&&'G'==tmp[i+1] ) {
            flag = 1;
            lengthNMEA = 0;
            //printf ( "\n\nflag -> 1, i = %d\n",i );
        }
        if ( ( u_char ) tmp[i] == 0xB5 && ( u_char ) tmp[i+1] == 0x62 ) {
            flag = 2;
            lengthUBX = 0;
            //printf ( "\n\nflag -> 2, i = %d\n",i );
        }
        if(1==flag){
            bufferNMEA[lengthNMEA] = tmp[i];
            lengthNMEA++;
        }
        if(2==flag){
            bufferUBX[lengthUBX] = tmp[i];
            if(6==lengthUBX){
                //length is writen in 4,5 in protocol
                lengthUBXProtocol = *(u_int16_t*)(bufferUBX+4);
            }
            lengthUBX++;
        }
        if(('\n'==tmp[i]||'\r'==tmp[i])&&1==flag){
            printf("Got a NMEA,l = %d.:",lengthNMEA);
//            if (showData)
            if (1)
                for(int k = 0;k<lengthNMEA;k++){
                    printf("%c",bufferNMEA[k]);
                }
            flag=0;
        }
        //there are 8 extra bytes besides the playload;
        if(lengthUBX==lengthUBXProtocol+8 && 2==flag){
//            printf("\nGot a UBX %02x %02x,l = %d, i = %d, count = %d\n",bufferUBX[2],
//                   bufferUBX[3], lengthUBXProtocol, i,count++);
            if(showData && lengthUBX<1024)
                for(int k = 0;k<lengthUBX;k++){
                    printf("%02x ",(u_char) bufferUBX[k]);
                }
            parse_UBX(bufferUBX);
            flag=0;
        }
    }
}

void SerialData::parse_UBX(char *buffer) {
    //todo:

    if(0x02==(u_char)buffer[2]){

        if(0x15==(u_char)buffer[3]){
            printf("\n--0--ParseRawData,%d\n",gnss->useBeiDou);
            gnss->ParseRawData(buffer, lengthUBX);
            printf("\n--1--ParseRawData,%d\n",gnss->useBeiDou);

        }
        if(0x13==(u_char)buffer[3]){
//            printf("\n----0----ParseBstSubFrame,%d\n",gnss->useBeiDou);
            gnss->svsManager.UpdateEphemeris(buffer);
//            printf("\n----1----ParseBstSubFrame,%d\n",gnss->useBeiDou);
        }
    }
}


