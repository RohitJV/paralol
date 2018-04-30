/*
Author : RohitJV
*/
#include <time.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>

#include "pr_graph.h"

#define TEST 1


// **************************************** Helper methods BEGIN ****************************************

/**
* @brief Return the number of seconds since an unspecified time (e.g., Unix
*        epoch). This is accomplished with a high-resolution monotonic timer,
*        suitable for performance timing.
*
* @return The number of seconds.
*/
static inline double monotonic_seconds() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/**
* @brief Output the seconds elapsed while execution.
*
* @param seconds Seconds spent on execution, excluding IO.
*/
static void print_time(double const seconds) {
  printf("Execution time: %0.04fs\n", seconds);
}
double startTime, endTime;

void print_graph(pr_graph * graph) {
	/*
	Print stuff to console.
	*/
	// printf("Vertices : %d\n", graph->nvtxs);
	// printf("Size of the object : %d\n", sizeof(*graph));
	
	/*
	Print stuff to file.
	*/
	FILE * opFile = fopen("test.txt", "a");		
	int i, start_vertex = graph->start_vertex;
	int edge_ptr = 0;
	for(i=1; i <= graph->nvtxs; i++) {		
		while (edge_ptr < graph->xadj[i]) {
			fprintf(opFile, "%d ", graph->nbrs[edge_ptr]+1);
			edge_ptr++;
		}
		fprintf(opFile, "\n");
	}	
	fclose(opFile);
}

void free_graph(pr_graph * graph) {
	free(graph->xadj);
	free(graph->nbrs);
	free(graph);
}

void send_graph(pr_graph * graph, int cur_proc_rank) {
	printf("Sending to processor %d\n", cur_proc_rank);
	int total_no_proc, rank_of_the_proc;
	MPI_Comm_size(MPI_COMM_WORLD, &total_no_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_of_the_proc);
  	  	  	
  	MPI_Send( (void*)&(graph->nvtxs), 1,  MPI_LONG_LONG_INT, cur_proc_rank, 0, MPI_COMM_WORLD);
  	MPI_Send( (void*)&(graph->nedges), 1,  MPI_LONG_LONG_INT, cur_proc_rank, 0, MPI_COMM_WORLD);
  	MPI_Send( (void*)&(graph->start_vertex), 1,  MPI_LONG_LONG_INT, cur_proc_rank, 0, MPI_COMM_WORLD);
  	printf("sent vertices : %d\n", graph->nvtxs);
  	MPI_Send( (void*)graph->xadj, graph->nvtxs + 1, MPI_LONG_LONG_INT, cur_proc_rank, 0, MPI_COMM_WORLD);  	
  	MPI_Send( (void*)graph->nbrs, graph->nedges, MPI_UNSIGNED_LONG, cur_proc_rank, 0, MPI_COMM_WORLD);
}

// **************************************** Helper methods END... ****************************************


void distributeVertices(FILE * fin, pr_int nvtxs, pr_int nedges, pr_graph **last_proc_graph) {	
    pr_int cur_proc_rank, edge_count=0;
	int total_no_proc, rank_of_the_proc;
	MPI_Comm_size(MPI_COMM_WORLD, &total_no_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_of_the_proc);
  	
	/* determine how many vertices should each processor have */
	pr_int vertices_per_process = nvtxs / total_no_proc;
	char * line = malloc(1024 * 1024);
  	size_t len = 0; 
  	pr_int v = 0, prev_v = 0, prev_e = 0;  	  	  	

	#if TEST
  		FILE * opFile = fopen("test.txt", "w");
		fprintf(opFile, "%d %d\n", nvtxs, nedges);
		fclose(opFile);
  	#endif

  	/* Divide the vertices */
  	for(cur_proc_rank = 0; cur_proc_rank < total_no_proc; cur_proc_rank++) {
  		/*
  		Set the flag to true; if it true after the while loop, then that means the graph hasn't yet been reallocated
  		*/
  		int flag = 1;

  		/* 
	  	Initialize the graph
	  	*/	  	
	  	pr_graph * graph = malloc(sizeof(*graph));
	  	graph->nvtxs = v;
	  	graph->nedges = nedges;
	  	graph->start_vertex = v;
	  	// TODO - Reallocate with smaller size
	  	graph->xadj = malloc((nvtxs + 1) * sizeof(*graph->xadj));
	  	graph->nbrs = malloc(nedges * sizeof(*graph->nbrs));
	  	
  		while(v < nvtxs) {
  			ssize_t read = getline(&line, &len, fin);		    
		    if(read == -1) {
		      fprintf(stderr, "ERROR: premature EOF at line %lu\n", v+1);		      
		      return;
		    }
			/* Check for sinks -- these make pagerank more difficult. */
		    if(read == 1) {
		      fprintf(stderr, "WARNING: vertex '%lu' is a sink vertex.\n", v+1);
		      continue;
		    }

			graph->xadj[v - prev_v] = edge_count - prev_e;

		    /* Foreach edge in line. */
		    char * ptr = strtok(line, " ");
		    while(ptr != NULL) {
		      char * end = NULL;
		      pr_int const e_id = strtoull(ptr, &end, 10);
		      /* end of line */
		      if(ptr == end)
		        break;		    
		      graph->nbrs[edge_count - prev_e] = e_id - 1; /* 1 indexed */ 	      
		      edge_count++;
		      ptr = strtok(NULL, " ");
		    }		    
		    v++;

		    /* 
		    If we have enough edges for a processor, break the loop, and communicate the vertices/edges to the respective processor
		    */
		    if( (cur_proc_rank == total_no_proc-1 && v==nvtxs) || (cur_proc_rank < total_no_proc-1 &&  v%vertices_per_process == 0)) {		    	
		    	/*
		    	Realloc the graph appropriately
		    	*/
		    	graph->nvtxs = v - prev_v;
		    	graph->nedges = edge_count - prev_e;
		    	graph->xadj = realloc(graph->xadj, (graph->nvtxs + 1) * sizeof(*graph->xadj));
		    	graph->xadj[graph->nvtxs] = graph->nedges;
		    	graph->nbrs = realloc(graph->nbrs, graph->nedges * sizeof(*graph->nbrs));
		    	
				if(cur_proc_rank != total_no_proc - 1) {
					send_graph(graph, cur_proc_rank);
					free_graph(graph);
				}
				else {
					*last_proc_graph = graph;
				}
		    	prev_e = edge_count;
		    	prev_v = v;
		    	flag = 0;
		    	break;
		    }		    
  		}  		
  	}
}

int main(int argc, char *argv[]) {
	if(argc!=3) {
		printf("Incorrect number of parameters passed\n");
		return 0;
	}

	MPI_Init(&argc, &argv);
	int total_no_proc, rank_of_the_proc, i=0;
	MPI_Comm_size(MPI_COMM_WORLD, &total_no_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_of_the_proc);

	int nvtxs, nedges;
	pr_graph * graph = malloc(sizeof(*graph));

	/*
	If the processor is the highest rank, then read the input file
	*/	
	if(rank_of_the_proc == total_no_proc-1) {
		// Make sanity checks
		if(argc == 1) {
		    fprintf(stderr, "usage: %s <graph> [output file]\n", *argv);
		    return EXIT_FAILURE;
		}
	    char * ifname = argv[1];
	  	char * ofname = NULL;
	  	if(argc > 2) {
	    	ofname = argv[2];
	  	}
	  	
	  	/* read nvtxs and nedges */
	  	FILE * fin = fopen(ifname, "r");
		fscanf(fin, "%d", &nvtxs);
		fscanf(fin, "%d", &nedges);			
		fscanf(fin, "\n"); /* make sure we process the newline, too. */

		printf("vertices: %d, edges: %d\n", nvtxs, nedges);	

	  	distributeVertices(fin, nvtxs, nedges, &graph);

	  	fclose(fin);
	}
	else {
		/*
		For each receiving processor, allocate space for graph
		*/		
		graph->nvtxs = nvtxs;
	  	graph->nedges = nedges;
	  	graph->start_vertex = 0;
	  	// TODO - Reallocate with smaller size	  	

	  	MPI_Status recv_status;
		MPI_Recv( (void*)&(graph->nvtxs), 1, MPI_LONG_LONG_INT, total_no_proc - 1, 0, MPI_COMM_WORLD, &recv_status);
		MPI_Recv( (void*)&(graph->nedges), 1, MPI_LONG_LONG_INT, total_no_proc - 1, 0, MPI_COMM_WORLD, &recv_status);
		MPI_Recv( (void*)&(graph->start_vertex), 1, MPI_LONG_LONG_INT, total_no_proc - 1, 0, MPI_COMM_WORLD, &recv_status);
		graph->xadj = malloc((graph->nvtxs + 1) * sizeof(*graph->xadj));
	  	graph->nbrs = malloc(graph->nedges * sizeof(*graph->nbrs));
	  	printf("Received vertices : %d\n", graph->nvtxs);
		MPI_Recv( (void*)(graph->xadj), graph->nvtxs + 1, MPI_LONG_LONG_INT, total_no_proc - 1, 0, MPI_COMM_WORLD, &recv_status);		
		MPI_Recv( (void*)(graph->nbrs), graph->nedges, MPI_UNSIGNED_LONG, total_no_proc - 1, 0, MPI_COMM_WORLD, &recv_status);
	}

	#if TEST
		for(i=0; i<total_no_proc; i++) {
			if(rank_of_the_proc == i)
				print_graph(graph);
			MPI_Barrier(MPI_COMM_WORLD);
		}				
	#endif 

	MPI_Finalize();
	return 0;
}