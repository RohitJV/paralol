#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 600
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */


#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include "pr_graph.h"


double monotonic_seconds() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

void print_time(double const seconds) {
  printf("Execution time: %0.04fs\n", seconds);
}


double startTime1, startTime2;

double * pagerank(
    pr_graph const * const graph,
    double const damping,
    int const max_iterations)
{
  /* grab graph structures to save typing */
  pr_int const nvtxs = graph->nvtxs;
  pr_int const * const restrict xadj = graph->xadj;
  pr_int const * const restrict nbrs = graph->nbrs;

  /* Initialize pageranks to be a probability distribution. */
  double * PR = malloc(nvtxs * sizeof(*PR));
  for(pr_int v=0; v < nvtxs; ++v) {
    PR[v] = 1. / (double) nvtxs;
  }

  /* Probability of restart */
  double const restart = (1 - damping) / (double) nvtxs;


  /* Convergence tolerance. */
  double const tol = 1e-12;

  double * PR_accum = malloc(nvtxs * sizeof(*PR));
  for(int i=0; i < max_iterations; ++i) {

    for(pr_int v=0; v < nvtxs; ++v) {
      PR_accum[v] = 0.;
    }

    /* Each vertex pushes PR contribution to all outgoing links */
    for(pr_int v=0; v < nvtxs; ++v) {
      double const num_links = (double)(xadj[v+1] - xadj[v]);
      double const pushing_val = PR[v] / num_links;

      for(pr_int e=xadj[v]; e < xadj[v+1]; ++e) {
        PR_accum[nbrs[e]] += pushing_val;
      }
    }

    /* Finalize new PR values */
    double norm_changed = 0.;
    for(pr_int v=0; v < nvtxs; ++v) {
      double const old = PR[v];
      PR[v] = restart + (damping * PR_accum[v]);

      norm_changed += (PR[v] - old) * (PR[v] - old);
    }
    norm_changed = sqrt(norm_changed);

    if(i > 1 && norm_changed < tol) {
      break;
    }
  }

  free(PR_accum);  
  return PR;
}


int main(
    int argc,
    char * * argv)
{
  startTime1 = monotonic_seconds();
  if(argc == 1) {
    fprintf(stderr, "usage: %s <graph> [output file]\n", *argv);
    return EXIT_FAILURE;
  }

  char * ifname = argv[1];
  char * ofname = NULL;
  if(argc > 2) {
    ofname = argv[2];
  }

  pr_graph * graph = pr_graph_load(ifname);  
  if(!graph) {
    return EXIT_FAILURE;
  }

  startTime2 = monotonic_seconds();

  double * PR = pagerank(graph, 0.85, 100);

  double endTime = monotonic_seconds();
  printf("Global Execution Time : %f\n", endTime - startTime1);
  printf("Algorithm Execution Time : %f\n", endTime - startTime2);
  
  /* write pagerank values */
  if(ofname) {
    FILE * fout = fopen(ofname, "w");
    if(!fout) {
      fprintf(stderr, "ERROR: could not open '%s' for writing.\n", ofname);
      return EXIT_FAILURE;
    }
    for(pr_int v=0; v < graph->nvtxs; ++v) {
      fprintf(fout, "%0.3e\n", PR[v]);
    }
    fclose(fout);
  }

  free(PR);  



 


  return EXIT_SUCCESS;
}