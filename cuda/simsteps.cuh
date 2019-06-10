/**
 * @file simsteps.cuh  
 * @author Russell Posner (rposner@uchc.edu)
 * @brief Contains most kernel code for sim
 * @version 0.1
 * @date 2019-06-09
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include <curand_kernel.h>
#include "cudaDmy.cuh"
#include "util.cuh"
#include "randomnumbers.cuh"


__device__ void apply_c_mat_as_char(const signed char * cmat, 
const int * ic_vec, int * rna_vec){
    unsigned int localtid = threadIdx.x;
    int myRowSum = 0;
    //Only the NMESS+nMicro threads participate
    if (localtid < NDNAM){
        const signed char * myRow = cmat+localtid*NMESS;
        for(int i=0;i<NMESS;i++){
            myRowSum += (int)myRow[i]*ic_vec[i];
        }
        if (myRowSum > 0){
            rna_vec[localtid]= 1;
        } else{
            rna_vec[localtid] = 0;
        }
    }
}





__device__ void apply_c_noise(curandStateMRG32k3a_t * states,  int * rna_vec,
const float cnoise0, const float cnoise1){
    INITTIDLOCALTID
    curandStateMRG32k3a_t * myState = &states[tid];
    if(localtid < NDNAM){
    float tmp = curand_uniform(myState);
    if (tmp < cnoise0)
        rna_vec[localtid]=0;
    else if (tmp < (cnoise0+cnoise1))
        rna_vec[localtid]=1;
    }
    
}


__device__ void generate_activemicro_vec(const int nActiveMicro, curandStateMRG32k3a_t * state,
int * vec){
    for(int i=0;i<NMICMAX;i++){
        vec[i] = (i<nActiveMicro) ? 1 : 0;
    }
    int i,j,tmp;
    for (i=NMICMAX-1;i>0;i--){
        j= curand(state) % ( i+1);
        tmp = vec[j];
        vec[j]=vec[i];
        vec[i]=tmp;
    }
}






__device__ void generate_l_probs(const float * lmat,const int * mir_vec,
const float rho, float * lprobs,const int * activeMicro){
    unsigned int localtid = threadIdx.x;
    //lprobs is of length NMESS
    if (localtid < NMESS){
    float mySuccessProb = 1.0f;
    const float * my_l_row = lmat+localtid*NMICMAX;
    for (int i=0;i<NMICMAX;i++){
        mySuccessProb *= (1.0f -(my_l_row[i] *mir_vec[i]*activeMicro[i]));
    }
    lprobs[localtid]=mySuccessProb*rho;
    }
}

__device__ void apply_l_probs(curandStateMRG32k3a_t * states, const float * lprobs,
const int * rna_vec, int * new_t_vec){
    INITTIDLOCALTID
    curandStateMRG32k3a_t * myState = &states[tid];
    if(localtid<NMESS){
        float tmp = curand_uniform(myState);
        new_t_vec[localtid] = rna_vec[localtid]*(tmp < lprobs[localtid]);
    }
}

//Really basic step implementation -> do TxC, do noise, do probs of inhibition, and then
//apply those and sum
__device__ void perform_step_cchar(curandStateMRG32k3a_t * states,const signed char * cmatchar,
 const float * Lmat,
const float Cnoise0,const float Cnoise1,float * lprobs,const float rho, int * activeMicro,
const int * ic_vec, int * rna_vec, int * new_state,int * new_sum){
    apply_c_mat_as_char(cmatchar,ic_vec,rna_vec);
    __threadfence_block();
    apply_c_noise(states,rna_vec,Cnoise0,Cnoise1);
    __threadfence_block();
    
    generate_l_probs(Lmat,rna_vec+NMESS,rho,
    lprobs,activeMicro);
    __threadfence_block();
    apply_l_probs(states,lprobs,rna_vec,new_state);
    __threadfence_block();
    vec_sum(new_state,NMESS,new_sum);

}



__device__ __inline__ void sharedLoad(int * src, int * dest){
    if(threadIdx.x < NMESS){
        dest[threadIdx.x] = src[threadIdx.x];
    }
    __threadfence_block();
}

__device__ __inline__ void clearTransmatRow(int * tmRow){
    if(threadIdx.x < (NMESS+1)){
        tmRow[threadIdx.x] = 0;
    }
    __threadfence_block();
}

__device__ __inline__ void storeTransmatRow(int *src, int * dest){
    if (threadIdx.x < NMESS){
        dest[threadIdx.x]=src[threadIdx.x];
    }
 
}

//each thread load a row
__device__ __inline__ void loadCMatShared(int * Csrc, signed char * cdst){
    INITTIDLOCALTID
    
    if(localtid < NDNAM){
        int * myRow = Csrc + NMESS * localtid;
        signed char * myDestRow = cdst + NMESS * localtid;
        for(int i=0;i<NMESS;i++){
            myDestRow[i]=(signed char) myRow[i];
        }
    }
    __threadfence_block();
}

///
__global__ void step_on_ic_vec(curandStateMRG32k3a_t * states, OpsData * oData,
int * ic_vec_array,int * transmatarray,const int nMicroToSet){
    __shared__ float my_lprobs[NMESS];
    __shared__ int my_rna_vec[NDNAM];
    __shared__ int my_new_state[NMESS];
    __shared__ int my_ic_vec[NMESS];
    __shared__ int myTransmatRow[NMESS+1];
    __shared__ signed char locCmat[NDNAM*NMESS]; 

    __shared__ int nActiveMicro;
    __shared__ int activeMicro[NMICMAX];

    int mySum;
    int * icCurr = ic_vec_array;
    INITTIDLOCALTID
    
    OpsData * my_oData = oData + blockIdx.x;
    int * myTransmat = transmatarray + TRANSMATDATASIZE* blockIdx.x;
    int * currTransmatRow = myTransmat;
    //Only thread 0 does shuffle, as parallel shuffle is a tricky one
    if(localtid == 0){
        nActiveMicro = nMicroToSet;
        generate_activemicro_vec(nActiveMicro,&states[tid],
        activeMicro);
    }
    loadCMatShared(my_oData->C,locCmat);
    __syncthreads();

    mySum = 0;
    
    for(int nOnes = 0; nOnes<=NMESS;nOnes++){
        clearTransmatRow(myTransmatRow);
        for(int i=0;i<NVPO;i++)
        {   
            sharedLoad(icCurr,my_ic_vec);
            mySum=0;
            perform_step_cchar(states,locCmat,my_oData->L,d_sParms.Cnoise0,d_sParms.Cnoise1,my_lprobs,d_sParms.rho,activeMicro,my_ic_vec,my_rna_vec,
            my_new_state,&mySum);
            //Need not be atomic as transmat row is previously loaded
            if(localtid == 0)
            {
                myTransmatRow[mySum]+=1;
            }
            __syncthreads();
            icCurr += NMESS;
        }
 
        storeTransmatRow(myTransmatRow,currTransmatRow);
        currTransmatRow += (NSTATES);
       __syncthreads();
    }

}


//A single thread uses a single shuffle to turn a (1,1,1,...1,0,..0) vector
//into a pseudorandom IC vector 
__device__ void generate_ic_vec(const int nOnes, curandStateMRG32k3a_t * state,
int * vec){
    for(int i=0;i<NMESS;i++){
        vec[i] = (i<nOnes) ? 1 : 0;
    }
    int i,j,tmp;
    for (i=NMESS-1;i>0;i--){
        j= curand(state) % ( i+1);
        tmp = vec[j];
        vec[j]=vec[i];
        vec[i]=tmp;
    }
}

//Supports IC generation kernel
__device__ void generate_mult_ics(const int nOnes, const int numVecs, 
curandStateMRG32k3a_t * states, int * vecs){
    INITTIDLOCALTID
    
    curandStateMRG32k3a_t * myState = &states[tid];
    const int numPerThread = numVecs / blockDim.x;
    const int remainderThreads = numVecs % blockDim.x;
    int * myVec = vecs+(NMESS * numPerThread * localtid);
    
    for(int i=0;i<numPerThread;i++){
        generate_ic_vec(nOnes,myState,myVec);
        myVec += NMESS;
    }
    if(localtid < remainderThreads){
        myVec = vecs+(NMESS * numPerThread * blockDim.x)+ (localtid*NMESS);
        generate_ic_vec(nOnes,myState,myVec);
    }
   
}


__global__ void generate_mult_ics_kernel(curandStateMRG32k3a_t * states,int * icstart){
    // if (gridDim.x == 1){
    //     generate_mult_ics(nOnes,numVecs,states,res);
    // }else{
        //each block will do some work on each value nOne - this tells us where to start each time
        const int numPerBlock = NVPO / gridDim.x;
        const int remainderVecs = NVPO % gridDim.x;
        const int myNumToGen = numPerBlock + 1*(blockIdx.x < remainderVecs);
        const int remainderLength = remainderVecs*myNumToGen * NMESS;
        const int baseOffset = blockIdx.x * myNumToGen * NMESS + (blockIdx.x >= remainderVecs)*remainderLength;
        int * startPosForNOne = NULL;
        int * myCurrStartPos = NULL;
        for(int i=0;i<(NMESS+1);i++){
            startPosForNOne = icstart + i * NVPO * NMESS;
            myCurrStartPos = startPosForNOne + baseOffset;
            generate_mult_ics(i,myNumToGen,states,myCurrStartPos);
        }
    //}
}

__device__ void generate_c_matrix_local(curandStateMRG32k3a_t * states,
 int * const resMat){
    INITTIDLOCALTID
    curandStateMRG32k3a_t * myState = &states[tid];
    if (localtid < NDNAM){
        //int rowsum=0;
        int * myRow = resMat+(localtid*NMESS);
        //rowsum = 
        generate_c_array(myState,myRow,NMESS);
    }
}

__device__ void generate_l_matrix_local(curandStateMRG32k3a_t * states,
float * const resLmat,float * mattmp){
    INITTIDLOCALTID
    curandStateMRG32k3a_t * myState = &states[tid];
    if (localtid < NMICMAX){
        float * myCol = mattmp + (NMESS * localtid);
        generate_l_column(myState,myCol);

        for(int i=0;i<NMESS;i++){
            resLmat[i*NMICMAX+localtid]=myCol[i];
        }
    }
}

//Each block builds its own "local" matrices which are then saved in the
//OpsArray
__global__ void mat_gen_kernel(curandStateMRG32k3a_t * states,
OpsData * const opsArray){
    extern __shared__ float lmattmp[];
    unsigned int localtid = threadIdx.x;
    OpsData * myOps = blockIdx.x + opsArray;
    int * local_cmat = myOps->C;
    float * local_lmat = myOps->L;
    generate_c_matrix_local(states,local_cmat);
    generate_l_matrix_local(states,local_lmat,lmattmp);

}