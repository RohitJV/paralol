/*!
\file  main.c
\brief This file is the entry point for paragon's various components
 
\date   Started 11/27/09
\author George
\version\verbatim $Id: omp_main.c 9585 2011-03-18 16:51:51Z karypis $ \endverbatim
*/


#include "simdocs.h"


int no_threads = 8;

/*************************************************************************/
/*! This is the entry point for finding simlar patents */
/**************************************************************************/
int main(int argc, char *argv[])
{
  params_t params;
  int rc = EXIT_SUCCESS;

  cmdline_parse(&params, argc, argv);

  printf("********************************************************************************\n");
  printf("sd (%d.%d.%d) Copyright 2011, GK.\n", VER_MAJOR, VER_MINOR, VER_SUBMINOR);
  printf("  nnbrs=%d, minsim=%.2f\n",
      params.nnbrs, params.minsim);

  gk_clearwctimer(params.timer_global);
  gk_clearwctimer(params.timer_1);
  gk_clearwctimer(params.timer_2);
  gk_clearwctimer(params.timer_3);
  gk_clearwctimer(params.timer_4);

  gk_startwctimer(params.timer_global);

  ComputeNeighbors(&params);

  gk_stopwctimer(params.timer_global);

  printf("    wclock: %.2lfs\n", gk_getwctimer(params.timer_global));
  printf("    timer1: %.2lfs\n", gk_getwctimer(params.timer_1));
  printf("    timer2: %.2lfs\n", gk_getwctimer(params.timer_2));
  printf("    timer3: %.2lfs\n", gk_getwctimer(params.timer_3));
  printf("    timer4: %.2lfs\n", gk_getwctimer(params.timer_4));
  printf("********************************************************************************\n");

  exit(rc);
}


/*************************************************************************/
/*! Reads and computes the neighbors of each document */
/**************************************************************************/
void ComputeNeighbors(params_t *params)
{
  int i, j, k, nhits;
  gk_csr_t *mat;
  int32_t *marker;
  gk_fkv_t *cand;
  FILE *fpout;

  printf("Reading data for %s...\n", params->infstem);

  mat = gk_csr_Read(params->infstem, GK_CSR_FMT_CSR, 1, 0);

  printf("#docs: %d, #nnz: %d.\n", mat->nrows, mat->rowptr[mat->nrows]);

  /* compact the column-space of the matrices */
  gk_csr_CompactColumns(mat);

  /* perform auxiliary pre-computations based on similarity */
  gk_csr_ComputeSquaredNorms(mat, GK_CSR_ROW);

  /* create the inverted index */
  gk_csr_CreateIndex(mat, GK_CSR_COL);

  /* create the output file */
  fpout = (params->outfile ? gk_fopen(params->outfile, "w", "ComputeNeighbors: fpout") : NULL);    

  /* find the best neighbors for each query document */
  gk_startwctimer(params->timer_1);

  // automate this part
  int div_x=2, div_y=4;  

  gk_fkv_t **total_hit_array = (gk_fkv_t **)malloc(mat->nrows * sizeof(gk_fkv_t*));
  total_hit_array[0] = (gk_fkv_t *)malloc(sizeof(gk_fkv_t) * mat->nrows * no_threads * params->nnbrs);
  for(i = 0; i < mat->nrows; i++)
    total_hit_array[i] = (*total_hit_array + no_threads * params->nnbrs * i);

  #pragma omp parallel for default(shared) private(i,j) num_threads(no_threads)
  for (i=0;i<no_threads;i++) {
    int offset1 = i/div_y, offset2 = i%div_y;
    int startRow1 = (mat->nrows)/div_x * offset1;
    int startRow2 = (mat->nrows)/div_y * offset2;
    int no_rows_y = (mat->nrows)/div_y;
    if(offset2 == div_y-1)
      no_rows_y = mat->nrows - startRow2;
    gk_csr_t *comparing_rows_mat = gk_csr_ExtractSubmatrix(mat, startRow2, no_rows_y);
    gk_csr_CreateIndex(comparing_rows_mat, GK_CSR_COL);    

    int endRow1 = startRow1 + (mat->nrows)/div_x;
    if(offset1 == div_x-1)
      endRow1 = mat->nrows;    

    for(j=startRow1;j<endRow1;j++) {      
      gk_fkv_t *hits = gk_fkvmalloc(comparing_rows_mat->nrows, "ComputeNeighbors: hits");
      int localhits = gk_csr_GetSimilarRows(comparing_rows_mat, 
                   mat->rowptr[j+1] - mat->rowptr[j], 
                   mat->rowind + mat->rowptr[j], 
                   mat->rowval + mat->rowptr[j], 
                   GK_CSR_JAC, params->nnbrs, params->minsim, hits, 
                   NULL, NULL);          
      int proc_offset = i*params->nnbrs;
      for (k=0; k<localhits; k++) {
        hits[k].val = hits[k].val + startRow2;
        total_hit_array[j][proc_offset + k] = hits[k];
      }
    }
  }  
  
  for(i=0;i<mat->nrows;i++) {
    int count = no_threads * params->nnbrs;    
    gk_fkvsortd(count, total_hit_array[i]);
    /* write the results in the file */
    if (fpout) {
      for (j=0; j<params->nnbrs; j++) 
        fprintf(fpout, "%8d %8zd %.3f\n", i, total_hit_array[i][j].val, total_hit_array[i][j].key);
    }
  }

  gk_stopwctimer(params->timer_1);

  /* cleanup and exit */
  if (fpout) gk_fclose(fpout);

  gk_csr_Free(&mat);

  return;
}