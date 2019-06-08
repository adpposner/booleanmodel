#ifndef CALL_GEN_SLN_H__
#define CALL_GEN_SLN_H__

#include "parameters.h"

//void solveSteps(char * outpth,const int nMMin, const int nMMax, const int nMStep,
//    C_Array_Init_Parms_flt cp, L_Array_Init_Parms_flt lp, Simulation_Params_flt sp);
//SysGenSln initializeSysGenParms(int nMi, int devId, char * outpth,
//C_Array_Init_Parms_flt cp, L_Array_Init_Parms_flt lp, Simulation_Params_flt sp);
#ifdef __cplusplus
extern "C"{
#endif
void setupSolutionSystem(const char * outpth,int devId,
SysGenSln *mySln, L_Array_Init_Parms_flt lp, C_Array_Init_Parms_flt cp,Simulation_Params_flt sp);
void closeSolutionSystem(SysGenSln * sln);
void doForParms(char * outpth,int nMicro, SysGenSln *sln);
#ifdef __cplusplus
}
#endif

#endif