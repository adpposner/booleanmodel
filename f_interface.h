/**
 * @file f_interface.h
 * @author Russell Posner (rposner@uchc.edu)
 * @brief Main interface to ``solution engine'' - i.e. Fortran subs in 
 * solved folder
 * @version 0.1
 * @date 2019-06-08
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef F_INTERFACE_RP_H__
#define F_INTERFACE_RP_H__

#include <stdlib.h>
#include <stdio.h>
#include <float.h>

//prevents intellisense from complaining (with some important caveats!!!!)
#ifdef __INTELLISENSE__
typedef double __float128;
#endif

//Basic model parameters. Must match those defined in solved/modelparms.f90.
//__float128 is non-standard!!!! ICC uses this designation for 128-bit float,
//but also known on some compilers as ``_Quad'' and on others (be very careful!!)
//as ``long double,'' which IS standards-compliant, but IS NOT a 128-bit float in
//nearly all instances. On Mac OS using icc, ``long double'' IS in fact a 16-byte quad-precision
//float
typedef struct modelparams
{
    int nMess;
    int nMicro;
    //TF connection params
    int tf_deg_low;
    int tf_deg_high;
    __float128 tf_p_nz;
    //miRNA ``defect'' - see manuscript
    __float128 defect;
    //Called ``pmr'' in paper, ratio of activating-to-inhibiting TFs
    __float128 tmpr;
    //miRNA connection params
    int m_deg_high;
    int m_deg_low;
    __float128 m_p_nz;
    //TxL efficiency
    __float128 rho;
    //TxC noise
    __float128 noise_pzero;
    __float128 noise_pone;
} modelparams;

//Global workspace vars. These include:
//nmirtargs - probabilities of each miRNA having N targets
//hypData - 3D hypergeometric tensor - used as a lookup table
//instead of using ``hyp2'' subroutine every time
//pKillData - 3D probability tensor for number of mRNAs ``killed''
//by miRNAs
//messnoise & micronoise - matrices for effect of TxC noise on
//mRNA and miRNA
//IMPORTANT - these data are MEANINGLESS on the C side, as there
//is no guarantee that __float128 is interpreted the same as
//REAL(KIND=16) - they are merely placeholders, and could just as
//easily be char * or int * or whatever (but ideally not void *)
typedef struct myGlobalVars
{
    __float128 *nmirtargs;
    __float128 *hypData;
    __float128 *pKillData;
    __float128 *messnoise;
    __float128 *micronoise;
    int allocd;
    int initd;
} myGlobalVars;

//Solution data vars. These include:
//condent - a MESSMAXSTATES-length double array used for cdtl entropy
//meanexp - mean TF expression (long-term)
//stddev - std of TF expression (long-term)
//entrate - entropy rate
//nRows & nCols - just set to MESSMAXSTATES for now
//allocPartial - if only condent has been alloc'd
//allocComplete - if the following extras have been alloc'd:
//total_res - a MESSMAXSTATES x MESSMAXSTATES vector w/ total transition mtx
//statdist - a MESSMAXSTATES vector w/ stationary dist'n
//eigval - there no matter what, but stores eigenvalue associated with stable dist'n (should be 1.0d0)
typedef struct solutionData
{
    double *condent;
    double meanexp;
    double stddev;
    double entrate;
    double eigval;
    double *statdist;
    double *total_res;
    int nRows;
    int nCols;
    int allocdpartial;
    int allocdcomplete;
} solutionData;


//these functions are for testing & extending methods
int initGlbls(modelparams *myParams, myGlobalVars varData);
int checkGlbls(modelparams *myParams, myGlobalVars varData);
int solvesystem(modelparams *myParams, myGlobalVars *glblVars, solutionData *slnData);


//core methods - allocation, iteration, and freeing - pretty self-explanatory.
//nextIteration will do initialization on first run
int allocateGlblsAndSln(myGlobalVars *vars, solutionData *slnData, int allocAll, FILE **fpp,const char * filename);
void nextIteration(modelparams myParams, myGlobalVars glblVars, solutionData *slnData, FILE *fp);
int freeSystem(myGlobalVars *vars, solutionData *slnData, FILE **fpp);



#endif