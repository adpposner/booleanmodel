#ifndef UTIL_H_RP__
#define UTIL_H_RP__

#include <float.h>


// #ifdef DBL_DECIMAL_DIG
// 	#define RP_DBL_DIGS	(DBL_DECIMAL_DIG)
// #else
// 	#ifdef DECIMAL_DIG
// 		#define RP_DBL_DIGS (DECIMAL_DIG)
// 	#else
// 		#define RP_DBL_DIGS (DBL_DIG+3)
// 	#endif
// #endif

// #define DBLFMT "%.*e"


/**
 * note choice of numbers here. Performance is affected by all of
 *  NMESS, NMICMAX, and NDNAM. The numbers chosen here (particulary NMICMAX)
 * have been chosen to keep the number of total items < 192 in general,but
 * performance appeared to be best only going up to 40.
 * NOTE: making NMESS <= 127 makes a HUGE difference
 */
#define NMESS   130
#define NMICMAX  40  
#define NDNAM    (NMESS + NMICMAX)
// #define TPB 256
#define WARP_SIZE_DEFAULT_RP    32
#define NSTATES (NMESS+1)
#define TRANSMATDATASIZE NSTATES * NSTATES
#define NMICROMIN 0
#define NMICROMAX 32
#define NMICROSTEP 10
#define NSTEPS 1+((NMICROMAX - NMICROMIN)/NMICROSTEP)

#define MESSMAX_SOLVED  151

//Really important number - specifies the number of IC vectors with k 1's and nMess+1 - k 0's
//For low/high k's, most will be repeated but that's okay
//Must be high to see less common events
#define NVPO    5000
//Also very important, number of networks to test for given parameter sets
#define NUM_NETS_PER_PARAMSET   1000

//const int Cmat_mem_size = NMESS*NDNAM*sizeof(int);
const int Lmat_mem_size = NMESS*NMICMAX*sizeof(float);



#endif