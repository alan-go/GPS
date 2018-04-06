//
// Created by zli on 18-1-6.
//

#ifndef DM100_UBLOXM8N_H
#define DM100_UBLOXM8N_H

#include <iostream>
#include <string>
#include <thread>
#include <mutex>
#include <queue>
#include <chrono>
#include <unistd.h>
#include <boost/asio.hpp>
#include "dm100.h"

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
    std::string m_serial_port;
    bool m_get_3d_position;

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

            char buffer[1024], GPS_Buffer[512];
            int length = 0;
            int ind = 0;
            int flag = 0,flagNew = 0;//0:null,1:NMEA,2:UBX;

            while ( !stopMessageCapture() ) {
                char tmp[128];
                char UBX_Buffer[512], NMEA_Buffer[128];
                auto transferred = sp.read_some ( boost::asio::buffer ( tmp ) );
                if ( transferred <= 0 ) {
                    usleep ( 100 );
                    printf ( "serial port return 0\n" );
                    continue;
                }

                int start = 0, end = 0;
                if ( length>950 ) {
                    length = 0;
                    flag = 0;
                    flagNew = 0;
                }
                memcpy ( buffer + length, tmp, transferred );
                length += transferred;

                printf ( "\ntransferred = %d,length = %d, flag = %d\n",transferred,length,flag );

                //int startUBX=0,endUBX = 0,startNMEA = 0,endNMEA = 0;
                bool moveHead = false;
                for ( int i = length-transferred; i<length; i++ ) {
                    bool getNew = false;

                    if ( '$'==buffer[i]&&'G'==buffer[i+1] ) {
                        getNew = true;
                        flagNew = 1;
                        printf ( "\nnew flag = 1, i = %d\n",i );
                    }
                    if ( ( u_char ) buffer[i] == 0xB5 && ( u_char ) buffer[i+1] == 0x62 ) {
                        getNew = true;
                        flagNew = 2;
                        printf ( "\nnew flag = 2, i = %d\n",i );

                    }
                    if ( 1==flag && ( '\n'==buffer[i]||'\r'==buffer[i] ) ) {
                        printf ( "\nget a r/n, i = %d,start = %d, oldflag = %d\n",i,start,flag );

                        end = i;
						int l = end-start;
                        if ( l<128 ) {
                            printf ( "\ncopy NMEA : %d -- %d, l = %d \n",start,i-1,l );
                            memcpy ( NMEA_Buffer,buffer+start,l );
                            flagNew = 0;
							for(int k = 0;k<l;k++){
								printf("%c",NMEA_Buffer[k]);
							}
                        }
                    }
                    if ( getNew ) {
                        printf ( "getNew: start = %d, i  = %d, oldflag = %d\n",start,i,flag );
                        moveHead = true;
                        end = i-1;
						int l = end-start;
                        if ( end>start&&flag!=0 ) {
                            if ( 2==flag && l<512 ) {
								printf ( "\ncopy ubx-%02x : %d -- %d, l=%d\n",(u_char) buffer[start+2],start,i-1,l );
                                memcpy ( UBX_Buffer,buffer+start,l );
								for(int k = 0;k<l;k++){
									printf("%02x ",(u_char) UBX_Buffer[k]);
								}
                            }
                        }
                        start = i;
                        getNew = false;
                    }

                    flag = flagNew;

                }
                if ( moveHead ) {
					printf("\nmove: %d -> 0, rest = %d",start,length-start);
                    memcpy ( buffer,buffer+start,1024-start );
                    length-=start;
                    start = 0;
                }


            }
            sp.close();
        } catch ( ... ) {
            printf ( "failed to open serial port\n" );
            return;
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
