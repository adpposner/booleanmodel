
/**
 * @file util.cu
 * @author Russell Posner (rposner@uchc.edu)
 * @brief Vector sum methods & one device method
 * @version 0.1
 * @date 2019-06-09
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "util.cuh"
#include "parameters.h"


//I prefer to keep my the C and CUDA separate
//this function containing CUDA code goes here
void initDeviceParms(const System_Parameters * sysparms){
    checkCuda(cudaMemcpyToSymbol(d_cParms,&sysparms->c_Parms,sizeof(d_cParms),0,cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpyToSymbol(d_lParms,&sysparms->l_Parms,sizeof(d_lParms),0,cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpyToSymbol(d_sParms,&sysparms->s_Parms,sizeof(d_sParms),0,cudaMemcpyHostToDevice));
}


/**
 * @brief Block reduce subroutines from
 * https://devblogs.nvidia.com/faster-parallel-reductions-kepler/
 * 
 */
__inline__  __device__
int warpReduceSum(int val) {
  for (int offset = WARP_SIZE_DEFAULT_RP /2; offset > 0; offset /= 2) 
    val += __shfl_down_sync(0xffffffff,val, offset);
  return val;
}

__inline__ __device__ 
int blockReduceSum(int val){
    static __shared__ int shared[32];
    int lane = threadIdx.x % WARP_SIZE_DEFAULT_RP;
    int wid = threadIdx.x / WARP_SIZE_DEFAULT_RP;
    val = warpReduceSum(val);

    if(lane ==0) shared[wid]=val;
    __syncthreads();

    val = (threadIdx.x < blockDim.x / WARP_SIZE_DEFAULT_RP) ? shared[lane] : 0;
    if (wid ==0)val = warpReduceSum(val);
    return val;
}


__device__ void vec_sum(const int * vIn,const int len, int * sum){
    int tmp=0;
    if (threadIdx.x < len){
        tmp = vIn[threadIdx.x];}else{tmp=0;}
    tmp = blockReduceSum(tmp);
    if(threadIdx.x == 0){
        *sum = tmp;
    }
    //adtomic add not strictly required
    // if ((threadIdx.x & (WARP_SIZE_DEFAULT_RP-1)) == 0){
    //     atomicAdd(sum,tmp);
    // }
}

