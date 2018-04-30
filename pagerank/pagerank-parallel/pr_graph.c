

/* ensure we have `getline()` */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "pr_graph.h"



pr_graph * pr_graph_load(
    char const * const ifname)
{
  FILE * fin = fopen(ifname, "r");
  if(!fin) {
    fprintf(stderr, "ERROR: could not open '%s' for reading.\n", ifname);
    return NULL;
  }

  pr_graph * graph = malloc(sizeof(*graph));

  /* read nvtxs and nedges */
  fscanf(fin, "%lu", &(graph->nvtxs));
  fscanf(fin, "%lu", &(graph->nedges));
  fscanf(fin, "\n"); /* make sure we process the newline, too. */

  graph->xadj = malloc((graph->nvtxs + 1) * sizeof(*graph->xadj));
  graph->nbrs = malloc(graph->nedges * sizeof(*graph->nbrs));

  /* How many edges we have read. */
  pr_int edge_ptr = 0;

  char * line = malloc(1024 * 1024);
  size_t len = 0;

  /* Read in graph one vertex at a time. */
  for(pr_int v=0; v < graph->nvtxs; ++v) {
    ssize_t read = getline(&line, &len, fin);
    if(read == -1) {
      fprintf(stderr, "ERROR: premature EOF at line %lu\n", v+1);
      pr_graph_free(graph);
      return NULL;
    }

    /* Store the beginning of the adjacency list. */
    graph->xadj[v] = edge_ptr;

    /* Check for sinks -- these make pagerank more difficult. */
    if(read == 1) {
      fprintf(stderr, "WARNING: vertex '%lu' is a sink vertex.\n", v+1);
      continue;
    }

    /* Foreach edge in line. */
    char * ptr = strtok(line, " ");
    while(ptr != NULL) {
      char * end = NULL;
      pr_int const e_id = strtoull(ptr, &end, 10);
      /* end of line */
      if(ptr == end) {
        break;
      }
      assert(e_id > 0 && e_id <= graph->nvtxs);

      graph->nbrs[edge_ptr++] = e_id - 1; /* 1 indexed */
      ptr = strtok(NULL, " ");
    }
  }
  assert(edge_ptr == graph->nedges);
  graph->xadj[graph->nvtxs] = graph->nedges;

  free(line);

  return graph;
}


void pr_graph_free(
    pr_graph * const graph)
{
  free(graph->xadj);
  free(graph->nbrs);
  free(graph);
}


