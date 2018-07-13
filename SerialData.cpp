//
// Created by root on 7/13/18.
//

#include "SerialData.h"

SerialData::SerialData() {

}

SerialData::~SerialData() {

}

void SerialData::StartCapture(const std::string serialPort, unsigned int baudRate) {
    while (!stopCapture){
        printf("capture\n");
        sleep(2);
    }
}
