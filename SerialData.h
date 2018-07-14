//
// Created by root on 7/13/18.
//

#ifndef SERIALDATA_H
#define SERIALDATA_H


#include <string>
#include <unistd.h>
#include <boost/asio.hpp>



class GNSS;
class SerialData {
public:
    std::string serialPort_;
    unsigned int baudRate;
    bool stopCapture;
    bool showData;
    GNSS* gnss;

public:
    SerialData();
    ~SerialData();
    void StartCapture(const std::string serialPort, unsigned int baudRate);

    void ScanSerialData(char *tmp,int transferred);

    void parse_UBX(char * buffer);

    void RecordData(const std::string fileName);

    void ReadFromFile(const std::string fileName);
private:
    int flag,count;
    int lengthNMEA, lengthUBX;
    uint16_t *lengthUBXProtocol;
    char bufferNMEA[128], bufferUBX[1024];
private:

};


#endif //SERIALDATA_H
