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
    int type,id;
    std::string serialPort_;
    boost::asio::serial_port *sp_;
    unsigned int baudRate;
    bool stopCapture{0};
    bool showData{0},logOpen{1},paraseDara{0},sendGGA{0};
    char saveName[128];
    pthread_t thread;
    std::deque<Solution,Eigen::aligned_allocator<Eigen::Vector3d>> solRaw;
    GNSS* gnss;

public:
    SerialData();
    SerialData(int type,int id,std::string portNane, unsigned int baudRate,bool logOpen,bool pData,GNSS* gnss)
                : type(type),id(id),serialPort_(portNane),baudRate(baudRate),logOpen(logOpen),paraseDara(pData),gnss(gnss){};
    ~SerialData();
    void StartCapture(const std::string serialPort, unsigned int baudRate, char *saveName);

    int StopDevice();

    void ScanSerialData(char *tmp,int transferred);

    void parse_UBX(char * buffer);

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
