#ifndef NTRIPRTK_H
#define NTRIPRTK_H

class GNSS;
class NtripRTK{
public:
    GNSS *gnss;
public:
    NtripRTK();

    ~NtripRTK();

    void UpdatePosition();
private:

private:
};

#endif