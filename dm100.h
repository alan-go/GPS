//
// Created by zli on 18-1-7.
//

// #include <opencv2/opencv.hpp>


#ifndef DM100_DM100_H
#define DM100_DM100_H

#define DM100_MSG_GPS        1
#define DM100_MSG_IMU        2
#define DM100_MSG_IMAGE      3

#define DM100_GRAVITY        9.81007

class DM100Message {
public:
    int message_type;
    long int timestamp;
};

class IMUMessage : public DM100Message {
public:
    float gyro[3], accl[3], compass[3];
};

class ImageMessage : public DM100Message {
public:
//     cv::Mat image;
    int frame_id;
};

class GPSMessage : public DM100Message {
public:
    int gps_msg_type;		//1 GPGGA, 2 GPGSA, 3 GPRMC

    //GPRMC and GPGGA
    bool gps_available;
    double longitude, latitude, altitude, waterplane;
    double utc_time;
    int num_satellites;

    //GPGSA
    int location_type;		//1 invalid, 2 2D, 3 3D
    double PDOP, HDOP, VDOP;
};

long int currentUSecsSinceEpoch();

bool stopUblox = 0;
bool stopMessageCapture(){
	return stopUblox;
};
// void enqueueMessage(std::shared_ptr<DM100Message> msg);

#endif //DM100_DM100_H
