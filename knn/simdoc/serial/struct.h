/*! 
\file 
\brief Data-structure definitions

This file contains various data structures used by paragon 

\date Started 11/27/09
\author George
\version $Id: struct.h 9585 2011-03-18 16:51:51Z karypis $
*/

#ifndef _SIMPAT_STRUCT_H_
#define _SIMPAT_STRUCT_H_



/*************************************************************************/
/*! This data structure stores the various variables that make up the 
 * overall state of the system. */
/*************************************************************************/
typedef struct {
  int nnbrs;                    /*!< The maximum number of nearest grants to output */
  float minsim;                 /*!< The minimum similarity to use for keeping neighbors */ 

  int verbosity;                /*!< The reporting verbosity level */

  char *infstem;                /*!< The filestem of the input file */
  char *outfile;                /*!< The filename where the output will be stored */

  /* timers */
  double timer_global;
  double timer_1;
  double timer_2;
  double timer_3;
  double timer_4;
} params_t;


#endif 
