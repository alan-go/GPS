

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
        auto parse_gps_message=[&] ( char *GPS_Buffer, GPSMessage &msg ) -> bool {
            //printf("%s\n", GPS_Buffer); return false;
            bool res = false;
            if ( memcmp ( GPS_Buffer, "$GPRMC,", 7 ) == 0 || memcmp ( GPS_Buffer, "$GNRMC,", 7 ) == 0 ) {
                //printf("%s\n", GPS_Buffer);
                res = parse_GPRMC ( GPS_Buffer, msg );
            } else if ( memcmp ( GPS_Buffer, "$GPGSA,", 7 ) == 0 || memcmp ( GPS_Buffer, "$GNGSA,", 7 ) == 0 ) {
                //printf("%s\n", GPS_Buffer);
                res = parse_GPGSA ( GPS_Buffer, msg );
            } else if ( memcmp ( GPS_Buffer, "$GPGGA,", 7 ) == 0 || memcmp ( GPS_Buffer, "$GNGGA,", 7 ) == 0 ) {
                //printf("%s\n", GPS_Buffer);
                res = parse_GPGGA ( GPS_Buffer, msg );
            }
            return res;
        };

        try {
            boost::asio::io_service ios;
            boost::asio::serial_port sp ( ios, serialPort );
            sp.set_option ( boost::asio::serial_port::baud_rate ( 115200 ) );
            printf ( "successfully opened port %s\n", serialPort.c_str() );

            //char buffer[1024], GPS_Buffer[512];
            int flag = 0;//0:null,1:NMEA,2:UBX;
			char bufferUBX[2048], bufferNMEA[128];
			int lengthUBX=0,lengthNMEA = 0;
			u_int16_t* lengthUBXProtocol;
            while ( !stopMessageCapture() ) {
                char tmp[20480];
                auto transferred = sp.read_some ( boost::asio::buffer ( tmp ) );
                if ( transferred <= 0 ) {
                    usleep ( 100 );
                    printf ( "serial port return 0\n" );
                    continue;
                }

//                printf ( "transferred = %d , flag = %d\n",transferred, flag );

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
							printf("MM");
//							printf("\nlengthUBXProtocol = %d\n",*lengthUBXProtocol);
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
						printf("get a UBX %02x %02x,l = %d, i = %d, count = %d\n",bufferUBX[2],bufferUBX[3], lengthUBX, i,count++);
						if(showData && lengthUBX<512)
//						for(int k = 0;k<lengthUBX;k++){
//							printf("%02x ",(u_char) bufferUBX[k]);
//						}
						parse_UBX(bufferUBX);
						flag=0;
					}

                }

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

    bool parse_GPRMC ( char *GPS_Buffer, GPSMessage &info ) {
        char usefulBuffer[256], *subString = strstr ( GPS_Buffer, "," ), *subStringNext;
        if ( subString == nullptr ) {
            return false;
        }
        info.gps_msg_type = GPS_MESSAGE_GPRMC;
        printf ( "Get a GPRMC\n" );

        for ( int i = 0; i < 6; i++ ) {
            subString++;    //skip a ','
            if ( ( subStringNext = strstr ( subString, "," ) ) != nullptr ) {
                if ( subStringNext == subString ) {  //the field is empty
                    continue;
                }
                switch ( i ) {
                case 0:                    //UTC time
                    break;
                case 1:                    //GPS status, located or not?
                    info.gps_available = ( subString[0] == 'A' );
                    break;
                case 2:                    //latitude
                    usefulBuffer[0] = subString[0];
                    usefulBuffer[1] = subString[1];
                    usefulBuffer[2] = 0;
                    info.latitude = atof ( usefulBuffer ) + atof ( subString + 2 ) / 60.0;
                    printf ( "RMC latitude,longtitude = %lf, ",info.latitude );
                    printf ( "%lf\n ",info.longitude );

                    break;
                case 3:                    //N_S
                    break;
                case 4:                    //longitude
                    usefulBuffer[0] = subString[0];
                    usefulBuffer[1] = subString[1];
                    usefulBuffer[2] = subString[2];
                    usefulBuffer[3] = 0;
                    info.longitude = atof ( usefulBuffer ) + atof ( subString + 3 ) / 60.0;
                    printf ( "%lf\n ",info.longitude );

                    break;
                case 5:                    //East or west
                    break;
                default:
                    break;
                }
                subString = subStringNext;
            }
        }
        return true;
    }

    bool parse_GPGSA ( char *GPS_Buffer, GPSMessage &info ) {
        char usefulBuffer[256], *subString = strstr ( GPS_Buffer, "," ), *subStringNext;
        if ( subString == nullptr ) {
            return false;
        }
        info.gps_msg_type = GPS_MESSAGE_GPGSA;
        printf ( "Get a GPGSA\n" );

        for ( int i = 0; i < 17; i++ ) {
            subString++;    //skip a ','
            if ( ( subStringNext = strstr ( subString, "," ) ) != nullptr ||
                    ( subStringNext = strstr ( subString, "*" ) ) != nullptr ) {
                if ( subStringNext == subString ) {  //the field is empty
                    continue;
                }
                if ( i == 0 ) {  //manual or auto
                } else if ( i == 1 ) {  //1 invalid, 2: 2D, 3: 3D
                    memcpy ( usefulBuffer, subString, subStringNext - subString );
                    usefulBuffer[subStringNext - subString] = 0;
                    info.location_type = atoi ( usefulBuffer );
                    if ( info.location_type == 3 && !m_get_3d_position ) {
                        m_get_3d_position = true;
                        printf ( "GPS: get 3d positon ready, You can GO!\n" );
                    }
                } else if ( i >= 2 && i <= 13 ) {  //#id of satellite
                } else if ( i == 14 ) {      //PDOP
                    memcpy ( usefulBuffer, subString, subStringNext - subString );
                    usefulBuffer[subStringNext - subString] = 0;
                    info.PDOP = atof ( usefulBuffer );
                } else if ( i == 15 ) {      //HDOP
                    memcpy ( usefulBuffer, subString, subStringNext - subString );
                    usefulBuffer[subStringNext - subString] = 0;
                    info.HDOP = atof ( usefulBuffer );
                } else if ( i == 16 ) {      //VDOP
                    memcpy ( usefulBuffer, subString, subStringNext - subString );
                    usefulBuffer[subStringNext - subString] = 0;
                    info.VDOP = atof ( usefulBuffer );
                }
                subString = subStringNext;
            }
        }
        return true;
    }

    bool parse_GPGGA ( char *GPS_Buffer, GPSMessage &info ) {
        char usefulBuffer[256], *subString = strstr ( GPS_Buffer, "," ), *subStringNext;
        if ( subString == nullptr ) {
            return false;
        }
        info.gps_msg_type = GPS_MESSAGE_GPGGA;
        printf ( "Get a GPGGA\n" );

        for ( int i = 0; i < 11; i++ ) {
            subString++;    //skip a ','
            if ( ( subStringNext = strstr ( subString, "," ) ) != nullptr ) {
                if ( subStringNext == subString ) {  //the field is empty
                    continue;
                }
                switch ( i ) {
                case 0:        //UTC time
                    memcpy ( usefulBuffer, subString, subStringNext - subString );
                    usefulBuffer[subStringNext - subString] = 0;
                    info.utc_time = atof ( usefulBuffer );
                    break;
                case 1:        //latitude
                    usefulBuffer[0] = subString[0];
                    usefulBuffer[1] = subString[1];
                    usefulBuffer[2] = 0;
                    info.latitude = atof ( usefulBuffer ) + atof ( subString + 2 ) / 60.0f;
                    printf ( "GGA latitude,longtitude = %lf, ",info.latitude );
                    break;
                case 2:        //north or south
                    break;
                case 3:        //longitude
                    usefulBuffer[0] = subString[0];
                    usefulBuffer[1] = subString[1];
                    usefulBuffer[2] = subString[2];
                    usefulBuffer[3] = 0;
                    info.longitude = atof ( usefulBuffer ) + atof ( subString + 3 ) / 60.0f;
                    printf ( "%lf\n ",info.longitude );

                    break;
                case 4:        //east or west
                    break;
                case 5:        //GPS status, located or not?
                    info.gps_available = ( subString[0] != '0' );
                    break;
                case 6:        //number of satellites
                    usefulBuffer[0] = subString[0];
                    usefulBuffer[1] = subString[1];
                    usefulBuffer[2] = 0;
                    info.num_satellites = atoi ( usefulBuffer );
                    break;
                case 7:        //HDOP
                    memcpy ( usefulBuffer, subString, subStringNext - subString );
                    usefulBuffer[subStringNext - subString] = 0;
                    info.HDOP = atof ( usefulBuffer );
                    break;
                case 8:        //altitude
                    memcpy ( usefulBuffer, subString, subStringNext - subString );
                    usefulBuffer[subStringNext - subString] = 0;
                    info.altitude = atof ( usefulBuffer );
                    break;
                case 9:        //M, meter
                    break;
                case 10:        //reference water plane
                    memcpy ( usefulBuffer, subString, subStringNext - subString );
                    usefulBuffer[subStringNext - subString] = 0;
                    info.waterplane = atof ( usefulBuffer );
                    break;
                default:
                    break;
                }
                subString = subStringNext;
            }
        }
        return true;
    }
};

#endif //DM100_UBLOXM8N_H
