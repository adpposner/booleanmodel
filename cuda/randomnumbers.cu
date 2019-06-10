/**
 * @file randomnumbers.cu
 * @author Russell Posner (rposner@uchc.edu)
 * @brief CURAND-based backend for PRNG
 * @version 0.1
 * @date 2019-06-09
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "randomnumbers.cuh"

//params stored in constant memory
__constant__ C_Array_Init_Parms d_cParms;
__constant__ L_Array_Init_Parms d_lParms;
__constant__ Simulation_Params d_sParms;
__constant__ System_Parameters d_sysParms;


extern "C" __device__  int generateBoolean(curandStateMRG32k3a_t * state, float p){
    float tmp = curand_uniform(state);
    return (tmp < p) ? 1 : 0;
}



__host__ curandStateMRG32k3a_t * randstates_alloc(dim3 dimGrid,dim3 dimBlock){
    curandStateMRG32k3a_t * res = NULL;
    checkCuda(cudaMalloc((void **)&res,dimGrid.x*dimBlock.x*sizeof(*res)));
    return res;
}

__host__  void randstates_dealloc(curandStateMRG32k3a_t * toFree){
    checkCuda(cudaFree(toFree));

}


__global__ void initRNG(curandStateMRG32k3a_t * randstates, const unsigned int seed){
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;

    curand_init(seed,tid,0,&randstates[tid]);
}




/// generate operators, can be slow
extern "C" __device__ inline int getBinomial(curandStateMRG32k3a_t * state, const int d_lo, const int d_hi, float p){
    int i;
    int n = d_hi - d_lo;
    int res=d_lo;
    // printf("%s:%d:%f\n",__FILE__,__LINE__,p);
    for(i=0;i<n;i++)
        res += generateBoolean(state,p);
    
    return res;
}




extern "C" __device__ inline int getUniformInt(curandStateMRG32k3a_t * state, const int vmaxnoninclusive){
    return (int)((float)vmaxnoninclusive * (curand_uniform(state)));
}

extern "C" __device__ inline void kshufflearrayInt(curandStateMRG32k3a_t * state, int * arr, const int len){
    int j;
    int tmp;
    int n;
    for(n=len;n>1;n--){
        j=getUniformInt(state,n);
        tmp = arr[n-1];
        arr[n-1]=arr[j];
        arr[j]=tmp;
    }
}

extern "C" __device__ inline void kshufflearrayFlt(curandStateMRG32k3a_t * state, float * arr, const int len){
    int j;
    float tmp;
    int n;
    for(n=len;n>1;n--){
        j=getUniformInt(state,n);
        tmp = arr[n-1];
        arr[n-1]=arr[j];
        arr[j]=tmp;
    }
}

extern "C" __device__ int generate_c_array(curandStateMRG32k3a_t * state,int * arr, int len){
    int numNZ,numP;
    numNZ = numP = 0;
    numNZ=getBinomial(state,d_cParms.d_lo,d_cParms.d_hi,d_cParms.pnz);
    numP=getBinomial(state,0,numNZ,d_cParms.pprob);
    int i;
    for(i=0;i<numP;i++){
        arr[i]=1;
    }

    for(;i<numNZ;i++){
        arr[i]=-1;
    }
    for(;i<len;i++){
        arr[i]=0;
    }
    kshufflearrayInt(state,arr,len);
    return (2*numP - numNZ);
}

extern "C" __device__ void generate_l_column(curandStateMRG32k3a_t * state, float * arr){
    int ntarg = getBinomial(state,d_lParms.d_lo,d_lParms.d_hi,d_lParms.ptarg);
    float targval = d_lParms.miRStrength / (float)ntarg;
    int i;
    for(i=0;i<ntarg;i++)
        arr[i]=targval;
    for(;i<NMESS;i++){
        arr[i]=0.0;
    }

    kshufflearrayFlt(state,arr,NMESS);
}