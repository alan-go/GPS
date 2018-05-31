

#ifndef UBLOX_SOLVER_H
#define UBLOX_SOLVER_H

class SvInfo
{
public:

};

class UbloxSolver
{
private:
    double tow;
    double prMes;
    double


public:
    void ParseRawData(const char* message);
    void ParseBstSubFrame(const char* message);
};

#endif //UBLOX_SOLVER_H
