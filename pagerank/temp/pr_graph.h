#ifndef PAGERANK_GRAPH_H
#define PAGERANK_GRAPH_H

#include <stdint.h>

typedef uint64_t pr_int;


/**
* @brief Structure representing a sparse, undirected, and unweighted  graph.
*/
typedef struct
{
  /** The number of vertices in the graph. */
  pr_int nvtxs;
  /** The number of edges in the graph. */
  pr_int nedges;

  /** The sparsity structure of the adjacency list. Vertex v has outgoing edges
   *  xadj[v] (inclusive) to xadj[v+1] (exclusive). */
  pr_int * xadj;

  /** The IDs associated with every outgoing edge. All values will be in the range
   *  [0, nvtxs). */
  pr_int * nbrs;
} pr_graph;



/**
* @brief Load a graph from an input file.
*
* @param ifname The name of the file to read.
*
* @return  A pointer to the graph, or NULL on error.
*/
pr_graph * pr_graph_load(
    char const * const ifname);


/**
* @brief Free all memory allocated by `pr_graph_load()`.
*
* @param graph The graph to free.
*/
void pr_graph_free(
    pr_graph * const graph);


#endif
