#ifndef NTRIPRTK_H
#define NTRIPRTK_H

#include <iostream>
#include <string>
#include <unistd.h>

#include <arpa/inet.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <boost/asio.hpp>



class GNSS;
class NtripRTK{
public:
    GNSS *gnss;
    bool stopRTK;
    std::string serverIP_;
    unsigned short port_;
    int sock_;

public:
    NtripRTK();

    ~NtripRTK();

    void UpdateGGA();
    int SentGGA(const char *bufferGGA, int length);
    bool NtripLogin(const std::string &rtk_protocol);
    void RecvThread();
private:

private:
};

#endif