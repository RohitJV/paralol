/*!
\file  
\brief This file contains function prototypes

\date Started 1/18/07
\author George
\version\verbatim $Id: proto.h 9628 2011-03-23 21:15:43Z karypis $ \endverbatim
*/

#ifndef _SIMPAT_PROTO_H_
#define _SIMPAT_PROTO_H_

#ifdef __cplusplus
extern "C"
{
#endif

/* cmdline.c */
void cmdline_parse(params_t *ctrl, int argc, char *argv[]);


/* main.c */
void ComputeNeighbors(params_t *params);

/* util.c */
void gk_csr_CompactColumns(gk_csr_t *mat);
void *gk_malloc(size_t, char *);
void *gk_realloc(void *oldptr, size_t nbytes, char *msg);
void gk_free(void **ptr1,...);
void gk_FreeMatrix(void ***r_matrix, size_t ndim1, size_t ndim2);
void gk_AllocMatrix(void ***r_matrix, size_t elmlen, size_t ndim1, size_t ndim2);
void gk_errexit(int signum, char *f_str,...);
void errexit(char *f_str,...);
int gk_dfkvkselect(size_t n, int topk, gk_fkv_t *cand);
void gk_fkvsortd(size_t n, gk_fkv_t *base);
uintmax_t gk_GetCurMemoryUsed();
uintmax_t gk_GetMaxMemoryUsed();
char *gk_strdup(char *orgstr);
FILE *gk_fopen(char *fname, char *mode, const char *msg);
void gk_fclose(FILE *fp);
int gk_fexists(char *fname);
void gk_getfilestats(char *fname, gk_idx_t *r_nlines, gk_idx_t *r_ntokens, 
        gk_idx_t *r_max_nlntokens, gk_idx_t *r_nbytes);
gk_idx_t gk_getline(char **lineptr, size_t *n, FILE *stream);
double gk_WClockSeconds(void);
uintmax_t gk_GetCurMemoryUsed();
uintmax_t gk_GetMaxMemoryUsed();
void gk_csr_Init(gk_csr_t *mat);
gk_csr_t *gk_csr_Create();
void gk_csr_CompactColumns(gk_csr_t *mat);
void gk_csr_Normalize(gk_csr_t *mat, int what, int norm);
void gk_csr_ComputeSquaredNorms(gk_csr_t *mat, int what);
gk_csr_t *gk_csr_Read(char *filename, int format, int readvals, int numbering);
void gk_csr_FreeContents(gk_csr_t *mat);
void gk_csr_Free(gk_csr_t **mat);
void gk_csr_CreateIndex(gk_csr_t *mat, int what);
gk_csr_t *gk_csr_ExtractSubmatrix(gk_csr_t *mat, int rstart, int nrows);
int gk_csr_GetSimilarRows(gk_csr_t *mat, int nqterms, int *qind, float *qval, 
        int simtype, int nsim, float minsim, gk_fkv_t *hits, int *i_marker,
        gk_fkv_t *i_cand);


#ifdef __cplusplus
}
#endif

#endif 
