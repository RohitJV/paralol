/*!
\file  main.c
\brief This file is the entry point for paragon's various components
 
\date   Started 11/27/09
\author George
\version\verbatim $Id: omp_main.c 9585 2011-03-18 16:51:51Z karypis $ \endverbatim
*/


#include "simdocs.h"


int no_threads = -1;

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
  printf("  nnbrs=%d, minsim=%.2f, minsim=%d\n",
      params.nnbrs, params.minsim, params.nthreads);

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
  no_threads = params->nthreads;  

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

  /* Set split-up of rows */
  int div_x, div_y;
  if(no_threads == 1)
    div_x = 1;
  else if(no_threads == 2)
    div_x = 2;
  else if(no_threads == 4)
    div_x = 2;
  else if(no_threads == 8)
    div_x = 2; 
  div_y = no_threads/div_x;

  /* Create a global array for hits - total_hit_array[row][hits across all processors/threads] */
  gk_fkv_t **total_hit_array = (gk_fkv_t **)malloc(mat->nrows * sizeof(gk_fkv_t*));
  total_hit_array[0] = (gk_fkv_t *)malloc(sizeof(gk_fkv_t) * mat->nrows * no_threads * params->nnbrs);
  for(i = 0; i < mat->nrows; i++)
    total_hit_array[i] = (*total_hit_array + no_threads * params->nnbrs * i);

  gk_fkv_t **result_array = (gk_fkv_t **)malloc(mat->nrows * sizeof(gk_fkv_t*));
  result_array[0] = (gk_fkv_t *)malloc(sizeof(gk_fkv_t) * mat->nrows * no_threads * params->nnbrs);
  for(i = 0; i < mat->nrows; i++)
    result_array[i] = (*result_array + no_threads * params->nnbrs * i);

  /* find the best neighbors for each query document */
  gk_startwctimer(params->timer_1);

  #pragma omp parallel for default(shared) private(i,j,k) num_threads(no_threads)
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
      int proc_offset = i*params->nnbrs;      
      int localhits = gk_csr_GetSimilarRows(comparing_rows_mat, 
                   mat->rowptr[j+1] - mat->rowptr[j], 
                   mat->rowind + mat->rowptr[j], 
                   mat->rowval + mat->rowptr[j], 
                   GK_CSR_JAC, params->nnbrs, params->minsim, total_hit_array[j]+proc_offset, 
                   NULL, NULL);                         
    }
  }  
  
  #pragma omp parallel for default(shared) private(i,j) num_threads(no_threads)
  for(i=0;i<mat->nrows;i++) {    
    int *idx = (int*)malloc(no_threads*sizeof(int));
    for(j=0;j<no_threads;j++)
      idx[j] = j*params->nnbrs;
    int neighbor;
    for(neighbor=0;neighbor<params->nnbrs;neighbor++) {       
      int maxIdx = -1;
      float maxSim = -1;
      for(j=0;j<no_threads;j++) {        
        if(total_hit_array[i][idx[j]].key > maxSim) {
          maxSim = total_hit_array[i][idx[j]].key;
          maxIdx = j;
        }
      }       
      result_array[i][neighbor] = total_hit_array[i][idx[maxIdx]];
      int offset = (mat->nrows)/div_y;      
      result_array[i][neighbor].val = result_array[i][neighbor].val + offset * maxIdx%div_y;
      idx[maxIdx]++;       
    }
  }

  /* write the results in the file */
  for(i=0;i<mat->nrows;i++) {
    if (fpout) {
      for (j=0; j<params->nnbrs && result_array[i][j].key >= params->minsim; j++) 
        fprintf(fpout, "%8d %8zd %.3f\n", i, result_array[i][j].val, result_array[i][j].key);
    }
  }

  gk_stopwctimer(params->timer_1);

  /* cleanup and exit */
  if (fpout) gk_fclose(fpout);

  gk_csr_Free(&mat);

  return;
}