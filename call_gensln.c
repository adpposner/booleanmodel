/**
 * @file call_gensln.c
 * @author Russell Posner (rposner@uchc.edu)
 * @brief This is a demo of the C interface to the general solution system
 * @version 0.1
 * @date 2019-06-08
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "f_interface.h"
#include "generatevars.h"

#define N_TOTAL_RUNS 10

///Demo function to run a couple times
int demo(void){
    int szb;
    modelparams myParams;
     myParams.nMess = 130;
 myParams.nMicro = 20;
 myParams.tf_deg_low =1 ;
 myParams.tf_deg_high = 14;
 myParams.tf_p_nz = 0.201;
 myParams.defect = 0.3;
 myParams.tmpr = 3.0;
 myParams.m_deg_high = 26;
 myParams.m_deg_low = 3;
 myParams.m_p_nz = 0.666;
 myParams.rho = 0.85;
 myParams.noise_pzero = 0.01;
myParams.noise_pone = 0.01;

myGlobalVars glblvars = {0};
    solutionData slnData = {0};
FILE * fp = NULL;
    
    int status;
    if (status = allocateGlblsAndSln(&glblvars,&slnData,0,&fp,"results.txt") ) {
        fprintf(stderr,"alloc failed return code %d\n",status);
    }
   // initGlbls(&myParams,glblvars);
    //checkGlbls(&myParams,glblvars);
    nextIteration(myParams,glblvars,&slnData,fp);
    myParams.nMicro = 10;
    nextIteration(myParams,glblvars,&slnData,fp);
    myParams.nMicro = 20;
    nextIteration(myParams,glblvars,&slnData,fp);

    freeSystem(&glblvars,&slnData,&fp);

}

//Main subroutine to regen - nTotalRuns is total runs 
//(including nMicro varying), fnameout is output filename
int runSolver(const int nTotalRuns,const char * fnameout){
    FILE * fpout = NULL;
    myGlobalVars glblVars = {0};
    solutionData slnData = {0};
    modelparams params = {0};
    int status,i;
    if (status = allocateGlblsAndSln(&glblVars,&slnData,0,&fpout,fnameout)){
        fprintf(stderr,"alloc failed return code %d\n",status);
        //May not work depending on precise status code, but also does no harm
        freeSystem(&glblVars,&slnData,&fpout);
        exit(-1);
    }
    //assuming mkl init works okay
    init_param_selector();
    
    for(i=0;i<nTotalRuns;i++){
        params = getNextParams();
        nextIteration(params,glblVars,&slnData,fpout);
        printf("currIteration = %d\n",i);
    }


    end_param_selector();
    freeSystem(&glblVars,&slnData,&fpout);
    return 0;
}


//insert favorite function here
int main(void){
    int status;
    //status = demo();
    status = runSolver(20,"res2.txt");
    return status;
}