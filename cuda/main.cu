/**
 * @file main.cu
 * @author Russell Posner (rposner@uchc.edu)
 * @brief Though officially "main," it pulls in a lot
 *  - e.g. simsteps.cuh which is pure source
 * @version 0.1
 * @date 2019-06-09
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include <stdlib.h>
#include "util.cuh"

#include "parameters.h"

#include "simsteps.cuh"

#define ENTRY_VALUE_NAME(name) #name

#define DEV_HOST_MALLOC(dev,nItems,casttype)                                        \
    casttype * toRet = NULL;                                                        \
    const int nBytes = nItems*sizeof(*toRet);                                       \
    fprintf(stderr,"%s : Allocating %d bytes\n",__func__,nBytes);                   \
    if(dev){                                                                        \
        checkCuda(cudaMalloc((void**)&toRet,nBytes));                               \
        checkCuda(cudaMemset(toRet,0,nBytes));                                      \
    }else{                                                                          \
         toRet = (casttype *) calloc(nItems,sizeof(*toRet));                        \
        if(!toRet){                                                                 \
            fprintf(stderr,"%s:%d - allocation type %s size %d failed",             \
            __FILE__,__LINE__,ENTRY_VALUE_NAME(casttype),nBytes);                   \
            exit(-5);                                                               \
        }                                                                           \
    }                                                                               \
    return toRet;                                                                   
                                                                    

int * allocICs(const int numVecs,int dev){
    const int numElts = numVecs*NMESS;
    DEV_HOST_MALLOC(dev,numElts,int)

}

OpsData * allocOpsData(dim3 dimGrid, int dev){
    const int numElts = dimGrid.x;
    DEV_HOST_MALLOC(dev,numElts,OpsData)

}

int * allocTransmats(dim3 dimGrid,int dev){
    const int numElts = dimGrid.x * TRANSMATDATASIZE;
    DEV_HOST_MALLOC(dev,numElts,int)

}



//writes out C and L matrices for all generated networks
void writeCMat(FILE * fp, int * Cmat){
    int i,j;
    int rowsum;
    for(i=0;i<NDNAM;i++){
        rowsum=0;
        for(j=0;j<NMESS-1;j++){
            rowsum+= *Cmat;
            fprintf(fp,"%d\t",*Cmat);
            Cmat++;
        }
        fprintf(fp,"%d\n",*Cmat);
        Cmat++;
    }
}

void writeLMat(FILE * fp, float * Lmat){
    int i,j;
    //float rowsums[NMESS]={0.0};
    //float colsums[NMICMAX] = {0.0};
    for(i=0;i<NMESS;i++){
        for(j=0;j<NMICMAX-1;j++){
            //rowsums[i]+=*Lmat;
            //colsums[j]+= *Lmat;
            fprintf(fp,"%f\t",*Lmat);
            Lmat++;
        }
        fprintf(fp,"%f\n",*Lmat);
        Lmat++;
    }
}


//generates the networks. generates dGrid.x networks and puts them in d_ops
void run_mat_gen_kernel(curandStateMRG32k3a_t * d_randStates, OpsData * d_ops, dim3 dGrid, dim3 dBlock){


    mat_gen_kernel KERNEL_ARG3(dGrid,dBlock,Lmat_mem_size) (d_randStates,d_ops);
    checkCuda(cudaDeviceSynchronize());
    #ifdef DEBUG_WRITE_MATRICES
    OpsData * h_opsData = NULL;
    size_t copyBytes = dGrid.x * sizeof(*d_opsData);
    h_opsData = allocOpsData(dGrid,0);
    checkCuda(cudaMemcpy(h_opsData,d_opsData,copyBytes,cudaMemcpyDeviceToHost));
    cudaDeviceSynchronize();
    
        char outbuf[245];
    for (int i=0;i<dGrid.x;i++){
        sprintf(outbuf,"LMatrixOut.%d.txt",i);
        FILE * fout = fopen(outbuf,"w");
        float * myLmat = &h_opsData[i].L[0];
        writeLMat(fout,myLmat);
    }
    for (int i=0;i<dGrid.x;i++){
        sprintf(outbuf,"CMatrixOut.%d.txt",i);
        FILE * fout = fopen(outbuf,"w");
        int * myCmat = &h_opsData[i].C[0];
        writeCMat(fout,myCmat);
    }
    free(h_opsData);
    #endif

}



void allocateSystem(curandStateMRG32k3a_t ** d_rand, int ** d_ics, int ** d_transmats, int ** h_transmats, 
        OpsData ** d_ops, dim3 dGrid, dim3 dBlock,int seed ){
    #ifdef DEBUG_WRITE_MATRICES
      cudaEvent_t start,stop;
    checkCuda(cudaEventCreate(&start));
    checkCuda(cudaEventCreate(&stop));
    checkCuda(cudaEventRecord(start));
    #endif
    *d_rand = randstates_alloc(dGrid,dBlock);
    *d_ics = allocICs(NVPO*(NMESS+1),1);
    *d_transmats = allocTransmats(dGrid,1);
    *h_transmats = allocTransmats(dGrid,0);
    *d_ops = allocOpsData(dGrid,1);

    initRNG KERNEL_ARG2(dGrid,dBlock) (*d_rand,seed);
    #ifdef DEBUG_WRITE_MATRICES
    checkCuda(cudaDeviceSynchronize());
    checkCuda(cudaEventRecord(stop));
    float millis = 0;
    checkCuda(cudaEventSynchronize(stop));
    checkCuda(cudaEventElapsedTime(&millis,start,stop));
    checkCuda(cudaEventDestroy(start));
    checkCuda(cudaEventDestroy(stop));
    printf("Allocations took %f s\n",millis/1000.0f);
    #endif
}

///
void initSystemParameters(const System_Parameters * myParams,dim3 dGrid, dim3 dBlock, curandStateMRG32k3a_t * const d_rand, 
 int * const d_transmats, int * const h_transmats, OpsData * const devOps) 
{
    #ifdef TIME_KERNELS_RP
    cudaEvent_t start,stop;
    checkCuda(cudaEventCreate(&start));
    checkCuda(cudaEventCreate(&stop));
    checkCuda(cudaEventRecord(start));
    #endif
    initDeviceParms(myParams);
    checkCuda(cudaMemset(d_transmats,0,dGrid.x*TRANSMATDATASIZE*sizeof(*d_transmats)));
    memset(h_transmats,0,sizeof(*h_transmats)*TRANSMATDATASIZE*dGrid.x);
    if (myParams->nMicro < NMICMAX){
    run_mat_gen_kernel(d_rand,devOps,dGrid,dBlock);
    }else{
        fprintf(stderr,"nMicro exceeds NMICMAX, skipping item %d - %d\n",myParams->nMicro,NMICMAX);
    }
    #ifdef TIME_KERNELS_RP
    checkCuda(cudaEventRecord(stop));
    float millis = 0;
    checkCuda(cudaEventSynchronize(stop));
    checkCuda(cudaEventElapsedTime(&millis,start,stop));
    checkCuda(cudaEventDestroy(start));
    checkCuda(cudaEventDestroy(stop));
    printf("setup System took %f s\n",millis/1000.0f);
    #endif

}

///
void freeSystem(curandStateMRG32k3a_t * randStates, OpsData * d_ops, int * d_ics, int * d_transmats,int * h_transmats){
    randstates_dealloc(randStates);
    checkCuda(cudaFree(d_ops));
    checkCuda(cudaFree(d_ics));
    checkCuda(cudaFree(d_transmats));
    free(h_transmats);
}



///
void step_kern(dim3 dGrid,dim3 dBlock,curandStateMRG32k3a_t * d_randStates,OpsData * d_ops, int * d_ics,
int * d_transmats,int myNMicro){
       #ifdef TIME_KERNELS_RP
   cudaEvent_t start,stop;
    checkCuda(cudaEventCreate(&start));
    checkCuda(cudaEventCreate(&stop));
    checkCuda(cudaEventRecord(start));
    #endif
    //regenerate ICs
    generate_mult_ics_kernel KERNEL_ARG2(dGrid,dBlock) (d_randStates,d_ics);
    checkCuda(cudaDeviceSynchronize());
    step_on_ic_vec KERNEL_ARG2(dGrid,dBlock) (d_randStates,d_ops,d_ics,d_transmats,myNMicro);
        #ifdef TIME_KERNELS_RP
    checkCuda(cudaEventRecord(stop));
    checkCuda(cudaEventSynchronize(stop));
    float millis =0;
    checkCuda(cudaEventElapsedTime(&millis,start,stop));
    printf("step took %f s \n",millis/1000.0f);
    #endif
}

extern "C" int loop(dim3 dGrid, dim3 dBlock, curandStateMRG32k3a_t * const d_rand, 
int * const d_ics, int * const d_transmats, int * const h_transmats, OpsData * const devOps){
     
    System_Parameters * pcurr = nextParameters();
    if (pcurr == NULL){
        fprintf(stderr,"Params are null, maybe EOF?\n");
        if(isendoffile())
            return 1;
        else
            return -1;//sp,
    }
    //This does regen networks when nMicro changes, could be made a bit faster by retaining networks
    //when only nmicro changes or only pzero/pone change
    initSystemParameters(pcurr,dGrid, dBlock, d_rand, d_transmats, h_transmats, devOps);
    step_kern(dGrid,dBlock,d_rand,devOps,d_ics,d_transmats,pcurr->nMicro);
    checkCuda(cudaPeekAtLastError());
    checkCuda(cudaDeviceSynchronize());
    checkCuda(cudaMemcpy(h_transmats,d_transmats,dGrid.x*TRANSMATDATASIZE*sizeof(*d_transmats),cudaMemcpyDeviceToHost));
    writeData(h_transmats,dGrid.x);
    return 0;
}






int main(int argc, const char * argv[]){
    int devId=0;
    //Number
    dim3 dGrid = {NUM_NETS_PER_PARAMSET,1,1};
    dim3 dBlock = {256,1,1};
    
    int numLoops;
    time_t t;
    srand((unsigned) time(&t));
    int myseed = rand();
    curandStateMRG32k3a_t * d_randStates;
    checkCuda( cudaSetDevice(devId) );

    int * d_ics, *d_transmats, *h_transmats;
    OpsData * device_Ops;

    if (argc < 4){
        fprintf(stderr,"Usage ./main inputfilename outputfileprefix numtrials");
        exit(-10);
    }else if ((numLoops = atoi(argv[3])) == 0){
        fprintf(stderr,"Invalid number of loops to run: %s\n",argv[3]);
        exit(-20);
    }
    int status;
    //setup
    printf("okay\n");
    openInputOutputFiles(argv[1],argv[2]);
    allocateSystem(&d_randStates,&d_ics,&d_transmats,&h_transmats,&device_Ops,dGrid,dBlock,myseed);
    //loop
    for(;numLoops > 0;numLoops--){
        status = loop(dGrid,dBlock,d_randStates,d_ics,d_transmats,h_transmats,device_Ops);
        if (status > 0) break;
        else if (status < 0) {
            fprintf(stderr,"Not eof, terminating\n");
            break;
        }
    }
    //shutdown
    freeSystem(d_randStates,device_Ops,d_ics,d_transmats,h_transmats);
    closeInputOutputFiles();
}
