/*!
\file   
\brief  This file contains the various header inclusions
 
\date   Started 11/27/09
\author George
\version\verbatim $Id: simdocs.h 9500 2011-03-03 15:42:05Z karypis $ \endverbatim
*/

#ifndef _SIMDOCS_HEADER_
#define _SIMDOCS_HEADER_

/*************************************************************************
* Header file inclusion section
**************************************************************************/
#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <signal.h>
#include <setjmp.h>
#include <assert.h>
#include <inttypes.h>
#include <sys/resource.h>
#include <sys/time.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


typedef ssize_t         gk_idx_t;         /* index variable */


/*-------------------------------------------------------------
 * The following data structure stores a sparse CSR format
 *-------------------------------------------------------------*/
typedef struct {
  int nrows, ncols;
  int *rowptr, *colptr, *rowids;
  int *rowind, *colind, *colids;
  float *rowval, *colval;
  float *rnorms, *cnorms;
  float *rsums, *csums;
} gk_csr_t;


#define SCNGKIDX "zd"

/* custom signals */
#define SIGMEM                  SIGABRT
#define SIGERR                  SIGABRT

/* CSR-related defines */
#define GK_CSR_ROW      1
#define GK_CSR_COL      2

#define GK_CSR_COS      1
#define GK_CSR_JAC      2
#define GK_CSR_MIN      3

#define GK_CSR_FMT_CSR          2
#define GK_CSR_FMT_CLUTO        1

#define LTERM                   (void **) 0     /* List terminator for GKfree() */



/*-------------------------------------------------------------
 * Program Assertions
 *-------------------------------------------------------------*/
#ifndef NDEBUG
#   define ASSERT(expr)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        assert(expr);                                                \
    }

#   define ASSERTP(expr,msg)                                          \
    if (!(expr)) {                                               \
        printf("***ASSERTION failed on line %d of file %s: " #expr "\n", \
              __LINE__, __FILE__);                               \
        printf msg ; \
        printf("\n"); \
        assert(expr);                                                \
    }
#else
#   define ASSERT(expr) ;
#   define ASSERTP(expr,msg) ;
#endif 


#define GK_MKKEYVALUE_T(NAME, KEYTYPE, VALTYPE) \
typedef struct {\
  KEYTYPE key;\
  VALTYPE val;\
} NAME;\

GK_MKKEYVALUE_T(gk_fkv_t,   float,    ssize_t);
GK_MKKEYVALUE_T(gk_dkv_t,   double,   ssize_t);


/*-------------------------------------------------------------
 * CSR conversion macros
 *-------------------------------------------------------------*/
#define MAKECSR(i, n, a) \
   do { \
     for (i=1; i<n; i++) a[i] += a[i-1]; \
     for (i=n; i>0; i--) a[i] = a[i-1]; \
     a[0] = 0; \
   } while(0) 

#define SHIFTCSR(i, n, a) \
   do { \
     for (i=n; i>0; i--) a[i] = a[i-1]; \
     a[0] = 0; \
   } while(0) 



#define gk_clearwctimer(tmr) (tmr = 0.0)
#define gk_startwctimer(tmr) (tmr -= gk_WClockSeconds())
#define gk_stopwctimer(tmr)  (tmr += gk_WClockSeconds())
#define gk_getwctimer(tmr)   (tmr)


#define gk_max(a, b) ((a) >= (b) ? (a) : (b))
#define gk_min(a, b) ((a) >= (b) ? (b) : (a))

#define QSSWAP(a, b, stmp) do { stmp = (a); (a) = (b); (b) = stmp; } while (0)



#include "defs.h"
#include "struct.h"
#include "proto.h"

#include "gk_getopt.h"
#include "gk_mksort.h"
#include "gk_mkmemory.h"
#include "gk_mkblas.h"



GK_MKALLOC_PROTO(gk_i,   int)
GK_MKALLOC_PROTO(gk_i32, int32_t)
GK_MKALLOC_PROTO(gk_f,   float)
GK_MKALLOC_PROTO(gk_fkv,   gk_fkv_t)
GK_MKALLOC_PROTO(gk_dkv,   gk_dkv_t)


#endif
