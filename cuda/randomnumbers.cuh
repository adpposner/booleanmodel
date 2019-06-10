/**
 * @file randomnumbers.cuh
 * @author Russell Posner (rposner@uchc.edu)
 * @brief Interface to generate pseudorandom numbers
 * @version 0.1
 * @date 2019-06-09
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef RANDOMNUMBERS_CUH_RP__
#define RANDOMNUMBERS_CUH_RP__
#include "util.cuh"
#include <curand_kernel.h>

#include "parameters.h"

extern "C"{

__device__  int generateBoolean(curandStateMRG32k3a_t * state, const float p);


__host__ curandStateMRG32k3a_t * randstates_alloc(dim3 dimGrid,dim3 dimBlock);

__host__  void randstates_dealloc(curandStateMRG32k3a_t * toFree);

__global__ void initRNG(curandStateMRG32k3a_t * randstates, const unsigned int seed);



__device__ int generate_c_array(curandStateMRG32k3a_t * state,int * arr, int len);
__device__ void generate_l_column(curandStateMRG32k3a_t * state, float * arr);

}

#endif