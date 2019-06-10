/**
 * @file parameters.h
 * @author Russell Posner (rposner@uchc.edu)
 * @brief Basic data structures with some I/O
 * @version 0.1
 * @date 2019-06-09
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef PARAMETERS_H_RP__
#define PARAMETERS_H_RP__

#include "util.h"



typedef struct C_Array_Init_Parms {
	int d_lo;
	int d_hi;
	float pnz;
	float pprob; // same as pmr / (1.0 + pmr)
	int byIndegree;
} C_Array_Init_Parms;

typedef struct L_Array_Init_Parms {
	int d_lo;
	int d_hi;
	float ptarg;
	float miRStrength; //same as 1 - defect
} L_Array_Init_Parms;

typedef struct Simulation_Params{
	float rho;
	float Cnoise0;
	float Cnoise1;
    //UNUSED
	float Lnoise0;
	float Lnoise1;
} Simulation_Params;

typedef struct System_Parameters{
    int nMess;
    int nMicro;
    C_Array_Init_Parms c_Parms;
    L_Array_Init_Parms l_Parms;
    Simulation_Params s_Parms;
} System_Parameters;

extern float accumulated[TRANSMATDATASIZE];



void openInputOutputFiles(const char * fnameIn,const char * fnameOut);
void closeInputOutputFiles();

extern System_Parameters * nextParameters();
extern int isendoffile();
extern void writeData(const int * h_tmats,int nMats);
//extern void initializeParametersXML(const char * cwd, System_Parameters * sys_params );


#endif