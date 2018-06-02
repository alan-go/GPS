#include "UbloxSolver.h"
using namespace std;


UbloxSolver::UbloxSolver(){}
UbloxSolver::~UbloxSolver(){}

void UbloxSolver::ParseRawData(const char *message) {
    cout<<"start solving rawdata"<<endl;
}

void UbloxSolver::ParseBstSubFrame(const char *message) {
    cout<<"Update subframe"<<endl;
}

