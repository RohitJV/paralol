/*!
\file  
\brief This file contains various macros
 
\date   Started 11/27/09
\author George
\version\verbatim $Id: macros.h 9496 2011-03-02 22:16:53Z karypis $ \endverbatim
*/

#ifndef _PARAGON_MACROS_H_
#define _PARAGON_MACROS_H_


/* Logfile printing macro */
#undef lprintf
#define lprintf(lvl, ...) \
  error_lprintf(lvl, __FILE__, __func__, __LINE__, __VA_ARGS__)

/* Message printing macro */
#undef mprintf
#define mprintf(lvl, ...) \
  error_mprintf(lvl, __FILE__, __func__, __LINE__, __VA_ARGS__)


/* A protected version of strcpy in which the size of the destination buffer
 * is checked for potential over-runs */
#define STRCPY(to, size, from)\
  do {\
    if (strlen(from) < size -1)\
      strcpy(to, from);\
    else {\
      strncpy(to, from, size-1);\
      to[size-1] = '\0';\
      lprintf(0, "strcpy overun: tosize: %zd, fromsize: %zd\n", size, strlen(from));\
    }\
  } while (0)




#define ACQUIRE_LOCK(lock)\
  do {\
    int rc;\
    if ((rc = pthread_mutex_lock(&(lock)))) \
      gk_errexit(SIGERR, "Failed on locking " # lock ": %s\n", gk_strerror(rc));\
  } while (0)

#define RELEASE_LOCK(lock)\
  do {\
    int rc;\
    if ((rc = pthread_mutex_unlock(&(lock)))) \
      gk_errexit(SIGERR, "Failed on unlocking " # lock ": %s\n", gk_strerror(rc));\
  } while (0)


#endif
  
 
