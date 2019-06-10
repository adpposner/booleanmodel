/**
 * @file util.cuh
 * @author Russell Posner (rposner@uchc.edu)
 * @brief Very basic **CUDA SPECIFIC** header, much like util.h
 * @version 0.1
 * @date 2019-06-09
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef UTIL_CUH_RP__
#define UTIL_CUH_RP__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cuda_runtime.h>

#include <device_launch_parameters.h>
#include "cudaDmy.cuh"
#include "util.h"
#include "parameters.h"



inline cudaError_t checkCuda_(cudaError_t result,const char filename[], int line){
    #if defined(RELEASE)
    return result;
    #else
    if (result != cudaSuccess){
        fprintf(stderr,"%s:%d CUDA Runtime Error: %s\n",filename,line,cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
    return result;
    #endif
}


#define checkCuda(res)  checkCuda_(res,__FILE__,__LINE__)

#define INITTIDLOCALTID unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x; unsigned int localtid = threadIdx.x;

__device__ void vec_sum(const int * vIn,const int len, int * sum);


//affects whether global values are used/not used - COMMON_CU must NOT be defined
// in this version
#ifndef COMMON_CU
extern __constant__ C_Array_Init_Parms d_cParms;
extern __constant__ L_Array_Init_Parms d_lParms;
extern __constant__ Simulation_Params d_sParms;
extern __constant__ System_Parameters d_sysParms;
#endif

typedef struct OpsData{
    int C[NDNAM*NMESS];
    float L[NMESS*NMICMAX];
} OpsData;


void initDeviceParms(const System_Parameters * sysparms);
__device__ int warpReduceSum(int val);


#endif