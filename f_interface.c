/*
f_interface.c - copyright Russell Posner 2018. 
This file contains all of the C subroutines needed to interface with the 
main solution contents (in the "solved" folder). This accompanies 
f_interface.h which contains all publicly exposed headers
 */

#include "f_interface.h"
#include <string.h>
//These values must be kept consistent between the Fortran and C sources
#ifndef MESSMAXSTATES
#define MESSMAXSTATES 151
#endif
#ifndef MICROMAXSTATES
#define MICROMAXSTATES 81
#endif
//basic debug print macro
#ifdef DEBUG_RP__
#define dbgprintf(...) fprintf(stderr, __VA_ARGS__)
#define WRITE_MATRICES 1
#else
#define dbgprintf(...) \
    do                 \
    {                  \
    } while (0)
#define WRITE_MATRICES 0
#endif

/////All of these functions are defined in solved/c_interface.f90///////////////

/*
Get sizes of workspace and solution vars in bytes so that they can be malloc'd on C side
 */
extern void get_sizes_bytes(int *mirtargprobsize,
                            int *hyptensorsize, int *pkilltensorsize, int *messnoisesize, int *micronoisesize,
                            int *condentsize, int *statdistsize, int *totalressize);

/*
Used to generate the workspace variables from myGlobalVars struct - note that the __float128 designations
here are meaningless on the C side, as there is poor quad support on most C compilers - it could just
as easily be any type
 */
extern void generate_global_vars(modelparams *myParams, __float128 *nmirtargsData,
                                 __float128 *pkillData, __float128 *hypsData, __float128 *noisemess, __float128 *noisemicro);
/*
This simply repeats generate_global_vars to make sure that the values persist properly. 
This is important - it does not test whether the values are ``correct,'' only that they persist 
 */
extern void check_global_vars(modelparams *myParams, __float128 *nmirtargsorig,
                              __float128 *pkillorig, __float128 *hypsorig);

/*
Main system solver routine - takes workspace and global vars, and returns 
double * condent - a [MESSMAXSTATES] length double vector with conditional entropies
double * meanexp -scalar value- a (ptr to) the mean active TF expression level
double * stddev -scalar value- the standard deviation of expression levels
double * entrate -scalar value- the entropy rate of the system
 */
extern void solve_sys_base(int *write_mats, modelparams *myParams, __float128 *nmirtargsData,
                           __float128 *pkills, __float128 *hyps, __float128 *noise_mess, __float128 *noise_micro, double *condent, double *meanexp,
                           double *stddev, double *entrate);

/*
Extended system solver routine - similar to above, but returns more info
double * condent - a [MESSMAXSTATES] vector with conditional entropies
double * statdist - a [MESSMAXSTATES] vector with the appropriate stable distribution
double * total_res - a [MESSMAXSTATES x MESSMAXSTATES] vector with the stochastic solution matrix
double * eigval - the eigenvalue associated wit hthe stationary dist'n
double * entrate - the entropy rate of the system
 */
extern void solve_sys_cplte(int *write_mats, modelparams *pset, __float128 *nmirtargsData,
                            __float128 *pkills, __float128 *hyps, __float128 *noise_mess, __float128 *noise_micro,
                            double *total_res, double *condent,
                            double *statdist, double *eigval, double *entrate);

/*
This routine takes the previous set of parameters used and the ``new'' parameters being used for the next 
iteration. It will compare the parameters, and regenerate the workspace data (myGlobalVars) that are needed
for the next iteration
 */
extern void update_global_vars(modelparams *params_old, modelparams *params_new, __float128 *nmirtargsData,
                               __float128 *pkillData, __float128 *hypsData, __float128 *noisemess, __float128 *noisemicro);

/////End of external (Fortran) function declarations///////////////

//Writes tab-delimited header to output file
void writeHeader(FILE *fp);

/*
allocateGlblsAndSln: responsible for allocation of workspace variables (myGlobalVars * vars) 
and for solution data (solutionData * slnData). 
Since there are two distinct solution subroutines (solve_sys_base and solve_sys_cplte), the 
Boolean value of allocAll tells whether to allocate the entire set of variables needed for sys_solve_cplte
or just the reduced form (used in manuscript). 
The FILE ** fpp is used to write the output data, and it is opened and header written to it
Returns a negative value if any allocs fail, 0 otherwise
 */
int allocateGlblsAndSln(myGlobalVars *vars, solutionData *slnData, int allocAll, FILE **fpp,const char * filename)
{
    __float128 *myHypData = NULL;
    __float128 *myPkillData = NULL;
    __float128 *myMirTargs = NULL;
    __float128 *myMessNoise = NULL;
    __float128 *myMicroNoise = NULL;
    double *mySysSln = NULL;
    double *myStatDist = NULL;
    double *myCondent = NULL;
    int hypsize, pkillsize, mirtargsize, condentsize, messnoisesize, micronoisesize, statdistsize, totalressize;

    hypsize = pkillsize = mirtargsize = condentsize = statdistsize = 0;
    messnoisesize = micronoisesize = totalressize = 0;
    get_sizes_bytes(&mirtargsize, &hypsize, &pkillsize, &messnoisesize, &micronoisesize, &condentsize, &statdistsize, &totalressize);
    dbgprintf("sizeof hyp = %d,sizeof pkill = %d,sizeof  mirtargs = %d, "
              "messnoise %d, micronoise %d\n",
              hypsize, pkillsize, mirtargsize, messnoisesize, micronoisesize);
    dbgprintf("sizeof condent = %d,sizeof statdist = %d,sizeof  total_res = %d\n", condentsize, statdistsize, totalressize);
    myHypData = malloc(hypsize);
    if (!myHypData)
    {
        fprintf(stderr, "ALLOC OF HYPDATA FAILED\n");
        return -1;
    }
    myPkillData = malloc(pkillsize);
    if (!myPkillData)
    {
        fprintf(stderr, "ALLOC OF PKILLDATA FAILED\n");
        free(myHypData);
        return -2;
    }
    myMirTargs = malloc(mirtargsize);
    if (!myMirTargs)
    {
        fprintf(stderr, "ALLOC OF MIRTARGS FAILED\n");
        free(myHypData);
        free(myPkillData);
        return -3;
    }
    myMessNoise = malloc(messnoisesize);
    if (!myMessNoise)
    {
        fprintf(stderr, "ALLOC OF MESSNOISE FAILED\n");
        free(myHypData);
        free(myPkillData);
        free(myMirTargs);
        return -5;
    }
    myMicroNoise = malloc(micronoisesize);
    if (!myMicroNoise)
    {
        fprintf(stderr, "ALLOC OF MiCRONOISE FAILED\n");
        free(myHypData);
        free(myPkillData);
        free(myMirTargs);
        free(myMessNoise);
        return -7;
    }

    myCondent = malloc(condentsize);
    if (!myCondent)
    {
        fprintf(stderr, "ALLOC OF CONDENT FAILED\n");
        free(myHypData);
        free(myPkillData);
        free(myMirTargs);
        free(myMessNoise);
        free(myMicroNoise);
        return -4;
    }

    if (allocAll)
    {
        myStatDist = malloc(statdistsize);
        if (!myStatDist)
        {
            fprintf(stderr, "ALLOC OF STATDIST FAILED\n");
            free(myHypData);
            free(myPkillData);
            free(myMirTargs);
            free(myCondent);
            free(myMessNoise);
            free(myMicroNoise);
            return -5;
        }
        mySysSln = malloc(totalressize);
        if (!mySysSln)
        {
            fprintf(stderr, "ALLOC OF SYSSLN FAILED\n");
            free(myHypData);
            free(myPkillData);
            free(myMirTargs);
            free(myCondent);
            free(mySysSln);
            free(myMessNoise);
            free(myMicroNoise);
            return -6;
        }
        slnData->allocdcomplete = 1;
    }
    slnData->allocdpartial = 1;
    dbgprintf("allocations completed\n");
    slnData->condent = myCondent;
    slnData->statdist = myStatDist;
    slnData->total_res = mySysSln;
    slnData->nRows = slnData->nCols = MESSMAXSTATES;
    vars->messnoise = myMessNoise;
    vars->micronoise = myMicroNoise;
    vars->hypData = myHypData;
    vars->pKillData = myPkillData;
    vars->nmirtargs = myMirTargs;
    vars->allocd = 1;
    vars->initd = 0;
    FILE *fp_temp = NULL;
    char fnamebuf[1000];
    memset(fnamebuf,0,1000);
    snprintf(fnamebuf,1000,"%s",filename);
    fp_temp = fopen(fnamebuf, "w");
    writeHeader(fp_temp);
    *fpp = fp_temp;
    return 0;
}

//Clears solution data for next iteration. Not typically completely necessary UNLESS nMess is changed!
//Run each iteration regardless
void clearSlnData(solutionData *slnData)
{
    memset(slnData->condent, 0, sizeof(*slnData->condent) * MESSMAXSTATES);
    slnData->eigval = slnData->entrate = slnData->meanexp = slnData->stddev = 0.0;

    if (slnData->allocdcomplete)
    {
        memset(slnData->statdist, 0, sizeof(*slnData->statdist) * MESSMAXSTATES);
        memset(slnData->total_res, 0, sizeof(*slnData->total_res) * MESSMAXSTATES * MESSMAXSTATES);
    }
}

//Can be used to set global vars for single run, but updateParameters works fine too
int initGlbls(modelparams *myParams, myGlobalVars varData)
{
    generate_global_vars(myParams, varData.nmirtargs, varData.pKillData, varData.hypData, varData.messnoise, varData.micronoise);
    return 0;
}

//Interface to the check_global_vars subroutine
int checkGlbls(modelparams *myParams, myGlobalVars varData)
{
    check_global_vars(myParams, varData.nmirtargs, varData.pKillData, varData.hypData);
    return 0;
}

//Used to compare parameters and change workspace variable values as needed.
//Just an interface to update_global_vars
int updateParameters(modelparams params_old, modelparams params_new, myGlobalVars varData)
{
    update_global_vars(&params_old, &params_new, varData.nmirtargs, varData.pKillData, varData.hypData,
                       varData.messnoise, varData.micronoise);
    return 0;
}

//Free global vars
int freeGlbls(myGlobalVars *vars)
{
    free(vars->hypData);
    free(vars->nmirtargs);
    free(vars->pKillData);
    free(vars->messnoise);
    free(vars->micronoise);
    vars->hypData = vars->nmirtargs = vars->pKillData = vars->messnoise = vars->micronoise = NULL;
    vars->allocd = vars->initd = 0;
}

//free solution data
int freeSlnData(solutionData *data)
{
    free(data->condent);
    free(data->statdist);
    free(data->total_res);
    data->condent = data->statdist = data->total_res = NULL;
    data->nCols = data->nRows = data->allocdcomplete = data->allocdpartial = 0;
    data->meanexp = data->stddev = data->entrate = data->eigval = 0.0;
    return 0;
}

void writeHeader(FILE *fp)
{
    fprintf(fp, "nMess\tnMicro\ttf_deg_low\ttf_deg_high\tm_deg_low\tm_deg_high\t");
    fprintf(fp, "tf_p_nz\tdefect\tpmr\tm_p_nz\trho\tp_zero\tp_one\t");
    fprintf(fp, "entropy_rate\texpression_mean\texpression_std");
    int i;
    for (i = 0; i < MESSMAXSTATES; i++)
        fprintf(fp, "\tcondent_%d", i);
    fputc('\n', fp);
}

//Writes output from the current iteration to fp - that includes params and slnData (base, not complete)
void writeIteration(FILE *fp, modelparams currParams, solutionData *slnData)
{
    fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t",
            currParams.nMess,
            currParams.nMicro,
            currParams.tf_deg_low,
            currParams.tf_deg_high,
            currParams.m_deg_low,
            currParams.m_deg_high);
    double tf_p_nz_tmp = currParams.tf_p_nz;
    double defect_tmp = currParams.defect;
    double tmpr_tmp = currParams.tmpr;
    double m_p_nz_tmp = currParams.m_p_nz;
    double rho_tmp = currParams.rho;
    double noise_pzero_tmp = currParams.noise_pzero;
    double noise_pone_tmp = currParams.noise_pone;
    fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t", tf_p_nz_tmp,
            defect_tmp,
            tmpr_tmp,
            m_p_nz_tmp,
            rho_tmp,
            noise_pzero_tmp,
            noise_pone_tmp);
    fprintf(fp, "%f\t%f\t%f\t", slnData->entrate, slnData->meanexp, slnData->stddev);
    int i;
    for (i = 0; i < MESSMAXSTATES - 1; i++)
    {
        fprintf(fp, "%f\t", slnData->condent[i]);
    }
    fprintf(fp, "%f\n", slnData->condent[i]);
}

//Main function exposed to iterate solutions.
//Takes params, workspace vars, slndata, and output FILE *
//will update the workspace vars, clear solution data, solve system and write output
void nextIteration(modelparams myParams, myGlobalVars glblVars, solutionData *slnData, FILE *fp)
{

    static modelparams orig_params = {0};
    //use as surrogate for undef'd
    if (orig_params.nMess == 0)
    {
        dbgprintf("First initialization\n");
        initGlbls(&myParams, glblVars);
        orig_params = myParams;
    }
    dbgprintf("updating...\n");
    updateParameters(orig_params, myParams, glblVars);
    clearSlnData(slnData);
    //checkGlbls(&myParams,glblVars);
    solvesystem(&myParams, &glblVars, slnData);
    dbgprintf("solved\n");
    writeIteration(fp, myParams, slnData);
    orig_params = myParams;
}

//Direct interface to system solver. Available in f_interface.h, but not intended to be used directly in
//most cases. Returns 1 if "cplte" solution runs successfully, 0 if "base" solution runs successfully,
//negative otherwise.
//Set
int solvesystem(modelparams *myParams, myGlobalVars *glblVars, solutionData *slnData)
{

    if (!slnData->allocdpartial)
    {
        fprintf(stderr, "%s:%d - Allocations not complete\n", __FILE__, __LINE__);
        return -1;
    }
    int write_mats = WRITE_MATRICES;
    if (slnData->allocdcomplete)
    {
        solve_sys_cplte(&write_mats, myParams, glblVars->nmirtargs, glblVars->pKillData, glblVars->hypData,
                        glblVars->messnoise, glblVars->micronoise,
                        slnData->total_res, slnData->condent, slnData->statdist, &slnData->eigval, &slnData->entrate);
        dbgprintf("completed, eigval = %f\tentrate = %f\n", slnData->eigval, slnData->entrate);
        return 1;
    }
    else
    {
        solve_sys_base(&write_mats, myParams, glblVars->nmirtargs, glblVars->pKillData, glblVars->hypData,
                       glblVars->messnoise, glblVars->micronoise, slnData->condent,
                       &slnData->meanexp, &slnData->stddev, &slnData->entrate);
        dbgprintf("completed, meanexp = %f\tstddev = %f\tentrate = %f\n", slnData->meanexp, slnData->stddev,
                  slnData->entrate);
        return 0;
    }
}

//Single subroutine used to free up vars and close solution file.
int freeSystem(myGlobalVars *vars, solutionData *slnData, FILE **fpp)
{
    freeGlbls(vars);
    freeSlnData(slnData);
    fclose(*fpp);
    *fpp = NULL;
    return 0;
}

#undef MESSMAXSTATES
#undef MICROMAXSTATES
#undef WRITE_MATRICES
#undef dbgprintf