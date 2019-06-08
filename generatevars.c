#include <mkl.h>
#include <mkl_vsl.h>
#include <math.h>
#include "generatevars.h"

#define NUM_TO_GEN 1000
#define N_DIM   11

//Be careful about bounds!!
#define N_MICRO_MIN 0
#define N_MICRO_MAX 79
#define N_MICRO_INTERVAL 13
#define N_MICRO_STEPS   1+((N_MICRO_MAX-N_MICRO_MIN) / N_MICRO_INTERVAL) 

#define TF_DEG_LOW_MIN  2
#define TF_DEG_LOW_MAX  12

#define TF_DEG_HIGH_MIN 6
#define TF_DEG_HIGH_MAX 40

#define TF_P_NZ_MIN 0.05
#define TF_P_NZ_MAX 0.40

#define DEFECT_MIN  0.1
#define DEFECT_MAX  0.3

#define TMPR_MIN    0.85
#define TMPR_MAX    4.10

#define M_DEG_HIGH_MIN  20
#define M_DEG_HIGH_MAX  70

#define M_DEG_LOW_MIN   3
#define M_DEG_LOW_MAX   20

#define M_P_NZ_MIN  0.1
#define M_P_NZ_MAX  0.6

#define RHO_MIN 0.40
#define RHO_MAX 0.99

#define NOISE_PZERO_MIN 0.01
#define NOISE_PZERO_MAX 0.25

#define NOISE_PONE_MIN 0.01
#define NOISE_PONE_MAX 0.25


static VSLStreamStatePtr stream;
static float rawParamsVector[NUM_TO_GEN][N_DIM];
static int ind = 0;
static int currNMicroStep = N_MICRO_MIN;
static modelparams currParams;

__float128 flt128inRange(const float prob,const float minval,const float maxval){
    __float128 toRet = (maxval - minval)*prob + minval;
    return toRet;
}

int intinRange(const float prob,const int minval,const int maxval){
    float tmp = (maxval - minval)*prob;
    int toRet = (int)round(tmp);
    toRet += minval;
    return toRet;
}

static modelparams flt_to_params(float * r11){
    modelparams toRet = {0};
    int tf_d_hi_tmp, mir_d_hi_tmp;
    toRet.nMess = 130;
    toRet.nMicro = 0;
    toRet.tf_deg_low = intinRange(r11[0],TF_DEG_LOW_MIN,TF_DEG_LOW_MAX);
    tf_d_hi_tmp = intinRange(r11[1],TF_DEG_HIGH_MIN,TF_DEG_HIGH_MAX);
    //make sure degrees are higher than low
    if (tf_d_hi_tmp <= toRet.tf_deg_low)
        toRet.tf_deg_high = toRet.tf_deg_low + (int)round(r11[1] * 10);
    else
        toRet.tf_deg_high = tf_d_hi_tmp;
    
    toRet.tf_p_nz =    flt128inRange(r11[2],TF_P_NZ_MIN,TF_P_NZ_MAX);
    toRet.defect = flt128inRange(r11[3],DEFECT_MIN,DEFECT_MAX);
    toRet.tmpr =   flt128inRange(r11[4],TMPR_MIN,TMPR_MAX);
    toRet.m_deg_low =  intinRange(r11[5],M_DEG_LOW_MIN,M_DEG_LOW_MAX);
    mir_d_hi_tmp = intinRange(r11[6],M_DEG_HIGH_MIN,M_DEG_HIGH_MAX);
    //make sure degrees are higher than low
    if (mir_d_hi_tmp <= toRet.m_deg_low)
        toRet.m_deg_high = toRet.m_deg_low + (int)round(10 * r11[6]);
    else
        toRet.m_deg_high = mir_d_hi_tmp;
    toRet.m_p_nz = flt128inRange(r11[7],M_P_NZ_MIN,M_P_NZ_MAX);
    toRet.rho =    flt128inRange(r11[8],RHO_MIN,RHO_MAX);
    toRet.noise_pzero =    flt128inRange(r11[9],NOISE_PZERO_MIN,NOISE_PZERO_MAX);
    toRet.noise_pone = flt128inRange(r11[10],NOISE_PONE_MIN,NOISE_PONE_MAX);
    return toRet;
}



void end_param_selector(){
    vslDeleteStream(&stream);
}

static modelparams getNext_(){
    if (ind == NUM_TO_GEN){
        vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,NUM_TO_GEN*N_DIM,
            (float*) rawParamsVector,0.0f,1.0f);
        ind = 0;
    }
    return flt_to_params(rawParamsVector[ind++]);
}

void init_param_selector(){
    vslNewStream(&stream,VSL_BRNG_SOBOL,N_DIM);
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,NUM_TO_GEN*N_DIM,
        (float*) rawParamsVector,0.0f,1.0f);
    ind=0;
    currParams = getNext_();
}


modelparams getNextParams(){
    if(currNMicroStep == N_MICRO_STEPS){
        currParams = getNext_();
        currNMicroStep = 0;
    }
    currParams.nMicro = currNMicroStep++ * N_MICRO_INTERVAL;
    printf("currparams nMess: %d\tnMicro: %d\n",currParams.nMess,currParams.nMicro);
    return currParams;
}

     
     
     
     
     
     
     
     
     
     