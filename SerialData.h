//
// Created by root on 7/13/18.
//

#ifndef SERIALDATA_H
#define SERIALDATA_H


#include <string>
#include <unistd.h>
#include <boost/asio.hpp>
#include "CommonInclude.h"



class GNSS;
class SerialData {
public:
    //mode:0:only care result; 1:raw data process; 2:
    int id,mode{0};
    //deviceType: 0:unknown,1:serialPort,2:file
    int deviceType{0};
    boost::asio::io_service ios;
    std::string serialPort_;
    boost::asio::serial_port *sp_;
    FILE* fpDev;
    unsigned int baudRate;
    bool stopCapture{0};
    bool showData{0},logOpen{1},sendGGA{0},ntripIn{0};
    char saveName[128];
    pthread_t thread;
    std::deque<Solution,Eigen::aligned_allocator<Eigen::Vector3d>> solRaw;
    GNSS* gnss;

public:
    SerialData();
    SerialData(int id,int type,std::string portNane, unsigned int baudRate,bool logOpen,GNSS* gnss)
            :id(id), mode(type),serialPort_(portNane),baudRate(baudRate),logOpen(logOpen),gnss(gnss){
        deviceType = 1;
    };
    SerialData(int id,int type,std::string fileNane,bool logOpen,GNSS* gnss)
            :id(id), mode(type),logOpen(logOpen),gnss(gnss){
        deviceType = 2;
        fpDev = fopen(fileNane.data(),"rb");
    };
    ~SerialData();
    void StartCapture(const std::string serialPort, unsigned int baudRate, char *saveName);

    int WtiteSerial(char* buffer, int len);

    int StopDevice();

    void ScanSerialData(char *tmp,int transferred);

    void parse_UBX(char * buffer);
    void parse_UBX_0106(char * buffer);

    int  ParaseGGA( char* gga);

    void RecordData(const std::string fileName);

    void ReadFromFile(const std::string fileName);
private:
    int flag,count;
    int lengthNMEA, lengthUBX;
    uint16_t lengthUBXProtocol;
    char bufferNMEA[128], bufferUBX[5120];
private:

};


#endif //SERIALDATA_H
