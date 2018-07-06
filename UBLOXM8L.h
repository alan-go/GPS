

#ifndef DM100_UBLOXM8N_H
#define DM100_UBLOXM8N_H

#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <mutex>
#include <queue>
#include <chrono>
#include <unistd.h>
#include <boost/asio.hpp>
#include "dm100.h"
#include "UbloxSolver.h"

#define GPS_MESSAGE_GPGGA       1
#define GPS_MESSAGE_GPRMC       2
#define GPS_MESSAGE_GPGSA       3


using namespace std;
class UBLOXM8L
{
public:
    UBLOXM8L() {
        m_thread = 0;
    }

    ~UBLOXM8L() {
    }

    int StartGPS ( const std::string &serialPort ) {
        m_serial_port = serialPort;
        m_get_3d_position = false;
        pthread_create ( &m_thread, nullptr, ThreadAdapter, this );
        return 1;
    }

    int StopGPS() {
        if ( m_thread ) {
            pthread_join ( m_thread, nullptr );
        }
        printf ( "quit GPS thread\n" );
        return 1;
    }

    void ScanSerialData(char *tmp,int transferred){
        int flag = 0;//0:null,1:NMEA,2:UBX;
        char bufferUBX[2048], bufferNMEA[128];
        int lengthUBX=0,lengthNMEA = 0;
        u_int16_t* lengthUBXProtocol;

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
                    lengthUBXProtocol = (u_int16_t*)(bufferUBX+4);
                    ////这里很奇怪，在mac上必须加一行打印，不然崩溃，Ubuntu不需要，还没找到原因。
//                  printf("MM");
//					printf("\nlengthUBXProtocol = %d\n",*lengthUBXProtocol);
                }
                lengthUBX++;
            }
            if(('\n'==tmp[i]||'\r'==tmp[i])&&1==flag){
                printf("get a NMEA,l = %d, i = %d\n",lengthNMEA,i);
                if (showData)
                    for(int k = 0;k<lengthNMEA;k++){
                        printf("%c",bufferNMEA[k]);
                    }
                flag=0;
            }
            //there are 8 extra bytes besides the playload;
            if(lengthUBX==*lengthUBXProtocol+8 && 2==flag){
                printf("\nget a UBX %02x %02x,l = %d, i = %d, count = %d\n",bufferUBX[2],bufferUBX[3], lengthUBX, i,count++);
                if(showData && lengthUBX<512)
                    for(int k = 0;k<lengthUBX;k++){
                        printf("%02x ",(u_char) bufferUBX[k]);
                    }
                parse_UBX(bufferUBX);
                flag=0;
            }

        }
    }

private:
    pthread_t m_thread;
    UbloxSolver solver;
    std::string m_serial_port;
    bool m_get_3d_position;
    bool showData  = true;
    int count = 0;

private:
    static void* ThreadAdapter ( void* __this ) {
        auto _this= ( UBLOXM8L* ) __this;
        _this->StartCapture ( _this->m_serial_port );
        return nullptr;
    }

    void StartCapture ( const std::string &serialPort ) {

        try {
            boost::asio::io_service ios;
            boost::asio::serial_port sp ( ios, serialPort );
            sp.set_option ( boost::asio::serial_port::baud_rate ( 115200 ) );
            printf ( "successfully opened port %s\n", serialPort.c_str() );

            while ( !stopMessageCapture() ) {
                char tmp[20480];
                auto transferred = sp.read_some ( boost::asio::buffer ( tmp ) );
                if ( transferred <= 0 ) {
                    usleep ( 100 );
                    printf ( "serial port return 0\n" );
                    continue;
                }
//                printf ( "transferred = %d , flag = %d\n",transferred, flag );
                ScanSerialData(tmp,transferred);
            }
            sp.close();
        } catch ( ... ) {
            printf ( "failed to open serial port\n" );
            return;
        }
    }

    void parse_UBX(char * buffer){
        if(0x02==(u_char)buffer[2]){
            if(0x15==(u_char)buffer[3]){
                printf("\n---ParseRawData\n");

//                solver.ParseRawData(buffer);
            }
            if(0x13==(u_char)buffer[3]){
                printf("\n-------ParseBstSubFrame\n");
                solver.ParseBstSubFrame(buffer);
            }
        }
    }


};

#endif //DM100_UBLOXM8N_H
