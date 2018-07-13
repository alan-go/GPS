//
// Created by root on 7/13/18.
//

#ifndef SERIALDATA_H
#define SERIALDATA_H


#include <string>
#include <unistd.h>


class GNSS;
class SerialData {
public:
    std::string serialPort_;
    unsigned int baudRate;
    bool stopCapture;
    GNSS* gnss;

public:
    SerialData();
    ~SerialData();
    void StartCapture(const std::string serialPort, unsigned int baudRate);
private:
private:

};


#endif //SERIALDATA_H
