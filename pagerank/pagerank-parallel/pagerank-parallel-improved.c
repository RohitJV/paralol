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
#include <assert.h>
#include "pr_graph.h"

#define TEST_CHUNKING 0
#define TEST_OUTGOING_TO_INCOMING 0
#define TEST_POST_COMMUNICATION 0

#define DEFAULT_LIST_SIZE 5000009


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

int cmpfunc (const void * a, const void * b) {
   return ( *(pr_int*)a > *(pr_int*)b ? 1 : -1);
}

void print_graph(pr_graph * graph) {
	/*
	Print stuff to console.
	*/
	// printf("Vertices : %d\n", graph->nvtxs);
	// printf("Size of the object : %d\n", sizeof(*graph));

	/*
	Print stuff to file.
	*/
	FILE * opFile = fopen("test_chunking.txt", "a");
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


void print_edges(pr_graph * graph) {
	/*
	Print stuff to console.
	*/
	// printf("Vertices : %d\n", graph->nvtxs);
	// printf("Size of the object : %d\n", sizeof(*graph));

	/*
	Print stuff to file.
	*/
	FILE * opFile = fopen("test_outgoing_to_incoming_main.txt", "a");
	int i, start_vertex = graph->start_vertex;
	int edge_ptr = 0;
	for(i=1; i <= graph->nvtxs; i++) {
		while (edge_ptr < graph->xadj[i]) {
			fprintf(opFile, "%d, %d\n", start_vertex + (i-1), graph->nbrs[edge_ptr]);
			edge_ptr++;
		}
	}
	fclose(opFile);
}

void free_graph(pr_graph * graph) {
	free(graph->xadj);
	free(graph->nbrs);
	free(graph);
}

void send_graph(pr_graph * graph, int cur_proc_rank) {
	int total_no_proc, rank_of_the_proc;
	MPI_Comm_size(MPI_COMM_WORLD, &total_no_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_of_the_proc);

  	MPI_Bcast( (void*)&(graph->total_nvtxs), 1, MPI_UNSIGNED, total_no_proc-1, MPI_COMM_WORLD );
  	MPI_Send( (void*)&(graph->nvtxs), 1,  MPI_UNSIGNED, cur_proc_rank, 0, MPI_COMM_WORLD);
  	MPI_Send( (void*)&(graph->nedges), 1,  MPI_UNSIGNED, cur_proc_rank, 0, MPI_COMM_WORLD);
  	MPI_Send( (void*)&(graph->start_vertex), 1,  MPI_UNSIGNED, cur_proc_rank, 0, MPI_COMM_WORLD);
  	MPI_Send( (void*)graph->xadj, graph->nvtxs + 1, MPI_UNSIGNED, cur_proc_rank, 0, MPI_COMM_WORLD);
  	MPI_Send( (void*)graph->nbrs, graph->nedges, MPI_UNSIGNED, cur_proc_rank, 0, MPI_COMM_WORLD);
}

pr_int calcOwnerProc(pr_int endpoint, pr_int total_nvtxs, int total_no_proc) {
	pr_int vertices_per_process = total_nvtxs / total_no_proc;
	if(endpoint/vertices_per_process < total_no_proc)
		return endpoint/vertices_per_process;
	return total_no_proc-1;
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

	#if TEST_CHUNKING
  		FILE * opFile = fopen("test_chunking.txt", "w");
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
	  	graph->total_nvtxs = nvtxs;
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

		    /* Foreach edge in line, add to a list */
		    pr_int* nnbrs_list = (pr_int *)malloc(DEFAULT_LIST_SIZE * sizeof(pr_int));
		    pr_int nnbrs_count = 0;
		    char * ptr = strtok(line, " ");
		    while(ptr != NULL) {
		      char * end = NULL;
		      pr_int const e_id = strtoull(ptr, &end, 10);
		      /* end of line */
		      if(ptr == end)
		        break;
		      nnbrs_list[nnbrs_count++] = e_id;
		      ptr = strtok(NULL, " ");
		    }
		    /* Sort this list and add to nbrs */
		    qsort(nnbrs_list, nnbrs_count, sizeof(pr_int), cmpfunc);
		    pr_int nnbrs_ptr = 0;
		    while(nnbrs_ptr != nnbrs_count) {
		      graph->nbrs[edge_count - prev_e] = nnbrs_list[nnbrs_ptr] - 1; /* 1 indexed */
		      edge_count++;
		      nnbrs_ptr++;
		    }
		    v++;
		    free(nnbrs_list);

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

  	free(line);
}

int main(int argc, char *argv[]) {
	if(argc!=3) {
		printf("Incorrect number of parameters passed\n");
		return 0;
	}

	double startTime = monotonic_seconds(), startTime2;

	MPI_Init(&argc, &argv);
	int total_no_proc, rank_of_the_proc, i=0, j=0, k=0;
	MPI_Comm_size(MPI_COMM_WORLD, &total_no_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_of_the_proc);

	int nvtxs, nedges, max_val = nvtxs+1000;
	pr_graph * graph = malloc(sizeof(*graph));
	char * ifname, * ofname;
	if(argc == 1) {
	    fprintf(stderr, "usage: %s <graph> [output file]\n", *argv);
	    return EXIT_FAILURE;
	}
    ifname = argv[1];
  	ofname = NULL;
  	if(argc > 2) {
    	ofname = argv[2];
  	}

	/*
	Read the input in chucks and distribute to processors
	*/
	if(rank_of_the_proc == total_no_proc-1) {
		/*
		If the processor is the highest rank, then read the input file
		*/

	  	/* read nvtxs and nedges */
	  	FILE * fin = fopen(ifname, "r");
		fscanf(fin, "%d", &nvtxs);
		fscanf(fin, "%d", &nedges);
		fscanf(fin, "\n"); /* make sure we process the newline, too. */

		// printf("vertices: %d, edges: %d\n", nvtxs, nedges);

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

	  	MPI_Status recv_status;

	  	MPI_Bcast( (void*)&(graph->total_nvtxs), 1, MPI_UNSIGNED, total_no_proc-1, MPI_COMM_WORLD );
		MPI_Recv( (void*)&(graph->nvtxs), 1, MPI_UNSIGNED, total_no_proc - 1, 0, MPI_COMM_WORLD, &recv_status);
		MPI_Recv( (void*)&(graph->nedges), 1, MPI_UNSIGNED, total_no_proc - 1, 0, MPI_COMM_WORLD, &recv_status);
		MPI_Recv( (void*)&(graph->start_vertex), 1, MPI_UNSIGNED, total_no_proc - 1, 0, MPI_COMM_WORLD, &recv_status);
		graph->xadj = malloc((graph->nvtxs + 1) * sizeof(*graph->xadj));
	  	graph->nbrs = malloc(graph->nedges * sizeof(*graph->nbrs));
		MPI_Recv( (void*)(graph->xadj), graph->nvtxs + 1, MPI_UNSIGNED, total_no_proc - 1, 0, MPI_COMM_WORLD, &recv_status);
		MPI_Recv( (void*)(graph->nbrs), graph->nedges, MPI_UNSIGNED, total_no_proc - 1, 0, MPI_COMM_WORLD, &recv_status);
	}
	MPI_Barrier(MPI_COMM_WORLD);


	/*
	Some random testing
	*/
	#if TEST_CHUNKING
		for(i=0; i<total_no_proc; i++) {
			if(rank_of_the_proc == i)
				print_graph(graph);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	#endif
	#if TEST_CHUNKING
		printf("%d\n", graph->total_nvtxs);
	#endif
	#if TEST_OUTGOING_TO_INCOMING || TEST_POST_COMMUNICATION
		if(rank_of_the_proc == total_no_proc-1) {
			FILE * opFile_test_main = fopen("test_outgoing_to_incoming_main.txt", "w");
			fprintf(opFile_test_main, ".......................................... \n");
			fclose(opFile_test_main);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		for(i=0; i<total_no_proc; i++) {
			if(rank_of_the_proc == i)
				print_edges(graph);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	#endif


	/*
	Setup to convert outgoing edges to incoming edges
	*/
	int* per_proc_count = (int*)calloc(total_no_proc, sizeof(int));
	int* per_proc_count_ptr = (int*)calloc(total_no_proc, sizeof(int));
	int* per_proc_count_start = (int*)calloc(total_no_proc, sizeof(int));
	int* per_proc_count_end = (int*)calloc(total_no_proc, sizeof(int));

	int* per_proc_outedges_count = (int*)calloc(total_no_proc, sizeof(int));
	int* per_proc_outedges_count_ptr = (int*)calloc(total_no_proc, sizeof(int));
	int* per_proc_outedges_count_start = (int*)calloc(total_no_proc, sizeof(int));
	int* per_proc_outedges_count_end = (int*)calloc(total_no_proc, sizeof(int));

	int edge_ptr=0, vtx_ptr;
	for(vtx_ptr = 1; vtx_ptr <= graph->nvtxs; vtx_ptr++) {
		int* per_proc_per_vtx_count = (int*)calloc(total_no_proc, sizeof(int));;
		while(edge_ptr < graph->xadj[vtx_ptr]) {
			pr_int endpoint = graph->nbrs[edge_ptr];
			int owner_proc_id = calcOwnerProc(endpoint, graph->total_nvtxs, total_no_proc);
			if(per_proc_per_vtx_count[owner_proc_id] == 0) {
				per_proc_per_vtx_count[owner_proc_id] = 1;
				per_proc_count[owner_proc_id]++;
			}
			per_proc_outedges_count[owner_proc_id]++;
			edge_ptr++;
		}
		free(per_proc_per_vtx_count);
	}

	int total_size = 0, grand_total_size=0;
	for(i=0; i<total_no_proc; i++) {
		per_proc_count_ptr[i] = total_size;
		per_proc_count_start[i] = total_size;
		total_size = total_size + per_proc_count[i];
		per_proc_count_end[i] = total_size;

		per_proc_outedges_count_ptr[i] = grand_total_size;
		per_proc_outedges_count_start[i] = grand_total_size;
		grand_total_size = grand_total_size + per_proc_outedges_count[i];
		per_proc_outedges_count_end[i] = grand_total_size;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	pr_int* source_nodes_for_proc = (pr_int*)malloc(total_size * sizeof(pr_int));
	pr_int* source_nodes_for_proc_count = (pr_int*)malloc(total_size  * sizeof(pr_int));
	pr_int* outedges_in_proc = (pr_int*)malloc(grand_total_size * sizeof(pr_int));
	edge_ptr = 0;
	for(vtx_ptr = 1; vtx_ptr <= graph->nvtxs; vtx_ptr++) {
		pr_int curr_actual_vtx = graph->start_vertex + (vtx_ptr-1);
		while(edge_ptr < graph->xadj[vtx_ptr]) {
			pr_int endpoint = graph->nbrs[edge_ptr];
			int owner_proc_id = calcOwnerProc(endpoint, graph->total_nvtxs, total_no_proc);
			int ptr = per_proc_count_ptr[owner_proc_id];
			if( ptr==per_proc_count_start[owner_proc_id]  ||  source_nodes_for_proc[ptr-1]!=curr_actual_vtx ) {
				source_nodes_for_proc_count[ptr] = 1;
				source_nodes_for_proc[ptr] = curr_actual_vtx;
				per_proc_count_ptr[owner_proc_id]++;
			}
			else {
				source_nodes_for_proc_count[ptr-1]++;
			}
			outedges_in_proc[ per_proc_outedges_count_ptr[owner_proc_id]++ ] = endpoint;
			edge_ptr++;
		}
	}

	free(per_proc_count_ptr);
	free(per_proc_count_start);
	free(per_proc_count_end);
	free(per_proc_outedges_count_ptr);
	free(per_proc_outedges_count_start);
	free(per_proc_outedges_count_end);

	/*
	Perform testing to check if the setup is consistent with the actual data
	*/
	#if TEST_OUTGOING_TO_INCOMING
		if(rank_of_the_proc == total_no_proc-1) {
			FILE * opFile_test = fopen("test_outgoing_to_incoming.txt", "w");
			fprintf(opFile_test, "....................................... \n");
			fclose(opFile_test);

		}
		MPI_Barrier(MPI_COMM_WORLD);

		for(i=0; i<total_no_proc; i++) {
			if(rank_of_the_proc == i) {
				FILE * opFile_test = fopen("test_outgoing_to_incoming.txt", "a");
				int vtx_ptr = 0, edg_ptr = 0;
				for(vtx_ptr=0; vtx_ptr<total_size; vtx_ptr++) {
					int temp_val = source_nodes_for_proc_count[vtx_ptr];
					while(temp_val>0) {
						fprintf(opFile_test, "%d, %d\n", source_nodes_for_proc[vtx_ptr], outedges_in_proc[edg_ptr]);
						temp_val--;
						edg_ptr++;
					}
				}
				fclose(opFile_test);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	#endif

	/**********************************************************************
		Perform exchange of information to convert outgoing to incoming
	***********************************************************************/
	int* per_proc_count_recv = (int*)calloc(total_no_proc, sizeof(int));
	int* per_proc_outedges_count_recv = (int*)calloc(total_no_proc, sizeof(int));
	MPI_Alltoall( (void*)per_proc_count, 1, MPI_UNSIGNED, (void*)per_proc_count_recv, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
	MPI_Alltoall( (void*)per_proc_outedges_count, 1, MPI_UNSIGNED, (void*)per_proc_outedges_count_recv, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

	int* per_proc_count_disp = (int*)calloc(total_no_proc, sizeof(int));
	int* per_proc_outedges_count_disp = (int*)calloc(total_no_proc, sizeof(int));
	int* per_proc_count_disp_recv = (int*)calloc(total_no_proc, sizeof(int));
	int* per_proc_outedges_count_disp_recv = (int*)calloc(total_no_proc, sizeof(int));
	/* Sending part */
	int startPosition = 0;
	for(i=0; i<total_no_proc; i++) {
		per_proc_count_disp[i] = startPosition;
		startPosition = startPosition + per_proc_count[i];
	}
	startPosition = 0;
	for(i=0; i<total_no_proc; i++) {
		per_proc_outedges_count_disp[i] = startPosition;
		startPosition = startPosition + per_proc_outedges_count[i];
	}
	/* Receiving part */
	startPosition = 0;
	for(i=0; i<total_no_proc; i++) {
		per_proc_count_disp_recv[i] = startPosition;
		startPosition = startPosition + per_proc_count_recv[i];
	}
	int total_size_recv = startPosition;
	startPosition = 0;
	for(i=0; i<total_no_proc; i++) {
		per_proc_outedges_count_disp_recv[i] = startPosition;
		startPosition = startPosition + per_proc_outedges_count_recv[i];
	}
	int grand_total_size_recv = startPosition;

	pr_int* source_nodes_for_proc_recv = (pr_int*)malloc(total_size_recv * sizeof(pr_int));
	MPI_Alltoallv( (void*)source_nodes_for_proc, (void*)per_proc_count, (void*)per_proc_count_disp, MPI_UNSIGNED, (void*)source_nodes_for_proc_recv, (void*)per_proc_count_recv, (void*)per_proc_count_disp_recv, MPI_UNSIGNED, MPI_COMM_WORLD);
	pr_int* source_nodes_for_proc_count_recv = (pr_int*)malloc(total_size_recv * sizeof(pr_int));
	MPI_Alltoallv( (void*)source_nodes_for_proc_count, (void*)per_proc_count, (void*)per_proc_count_disp, MPI_UNSIGNED, (void*)source_nodes_for_proc_count_recv, (void*)per_proc_count_recv, (void*)per_proc_count_disp_recv, MPI_UNSIGNED, MPI_COMM_WORLD);
	pr_int* outedges_in_proc_recv = (pr_int*)malloc(grand_total_size_recv * sizeof(pr_int));
	MPI_Alltoallv( (void*)outedges_in_proc, (void*)per_proc_outedges_count, (void*)per_proc_outedges_count_disp, MPI_UNSIGNED, (void*)outedges_in_proc_recv, (void*)per_proc_outedges_count_recv, (void*)per_proc_outedges_count_disp_recv, MPI_UNSIGNED, MPI_COMM_WORLD);


  for(i=0; i<grand_total_size_recv; i++) {
    outedges_in_proc_recv[i] -= graph->start_vertex;
  }

	#if TEST_POST_COMMUNICATION
		if(rank_of_the_proc == total_no_proc-1) {
			FILE * opFile_test = fopen("test_outgoing_to_incoming_post_comm.txt", "w");
			fprintf(opFile_test, "....................................... \n");
			fclose(opFile_test);

		}
		MPI_Barrier(MPI_COMM_WORLD);

		for(i=0; i<total_no_proc; i++) {
			if(rank_of_the_proc == i) {
				FILE * opFile_test = fopen("test_outgoing_to_incoming_post_comm.txt", "a");
				int vtx_ptr = 0, edg_ptr = 0;
				for(vtx_ptr=0; vtx_ptr<total_size_recv; vtx_ptr++) {
					int temp_val = source_nodes_for_proc_count_recv[vtx_ptr];
					while(temp_val>0) {
						fprintf(opFile_test, "%d, %d\n", source_nodes_for_proc_recv[vtx_ptr], outedges_in_proc_recv[edg_ptr]);
						temp_val--;
						edg_ptr++;
					}
				}
				fclose(opFile_test);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	#endif


	MPI_Barrier(MPI_COMM_WORLD);
	if(rank_of_the_proc==0) {
	  	double endTime = monotonic_seconds();
		printf("Execution Time till converting edges : %f\n", endTime - startTime);
		startTime2 = monotonic_seconds();
	}

	/**********************************************************************
		Page Rank !!!
	***********************************************************************/
	MPI_Request reqs[2];
	MPI_Status status;

	int max_iterations = 100;
	double const damping = 0.85;
	double *PR = malloc(graph->nvtxs * sizeof(double));
	pr_int v;
	for(v=0; v < graph->nvtxs; ++v)
    	PR[v] = 1. / (double) (graph->total_nvtxs);
	/* Probability of restart */
  	double const restart = (1 - damping) / (double) (graph->total_nvtxs);
  	/* Convergence tolerance. */
  	double const tol = 1e-9;

  	double* PR_accum;
  	double* PR_send = (double*)malloc(total_size * sizeof(double));
  	double* PR_recv = (double*)malloc(total_size_recv * sizeof(double));
  	for(i=0; i < max_iterations; ++i) {
  		if(rank_of_the_proc==0) {
  			printf("Iterations : %d\n", i);
  		}
  		pr_int v;
		PR_accum = calloc(sizeof(double), graph->nvtxs);

		// Replace the 'source_nodes_for_proc' elements with their pageranks
		for(j=0; j<total_size; j++) {
			pr_int local_vtx_idx = source_nodes_for_proc[j] - graph->start_vertex;
			PR_send[j] = PR[ local_vtx_idx ] / (double) (graph->xadj[local_vtx_idx+1] - graph->xadj[local_vtx_idx]);
		}

	 	/* Perform AllToAllv */
		MPI_Alltoallv( (void*)PR_send, (void*)per_proc_count, (void*)per_proc_count_disp, MPI_DOUBLE, (void*)PR_recv, (void*)per_proc_count_recv, (void*)per_proc_count_disp_recv, MPI_DOUBLE, MPI_COMM_WORLD);
		// for(j=1; j<=total_no_proc; j++) {
		// 	double s1,s2,e1,e2,s3,e3;
		// 	s1 = monotonic_seconds();
		// 	int send_pos = (rank_of_the_proc + j) % total_no_proc;
		// 	MPI_Isend(PR_send + per_proc_count_disp[send_pos], per_proc_count[send_pos], MPI_DOUBLE, send_pos, 1, MPI_COMM_WORLD, &reqs[0]);
		// 	int recv_pos = (rank_of_the_proc - j + total_no_proc) % total_no_proc;
		// 	MPI_Irecv(PR_recv + per_proc_count_disp_recv[recv_pos], per_proc_count_recv[recv_pos], MPI_DOUBLE, recv_pos, 1, MPI_COMM_WORLD, &reqs[1]);
    //
		// 	int cur_send_pos = (rank_of_the_proc + j - 1) % total_no_proc;
		// 	int cur_recv_pos = (rank_of_the_proc - j + 1 + total_no_proc) % total_no_proc;
    //
		// 	s2 = monotonic_seconds();
		// 	int k=0;
		// 	if(j==1) {
		// 		int edg_ptr = per_proc_outedges_count_disp[cur_send_pos];
		// 		for(k=0; k<per_proc_count[cur_send_pos]; k++) {
		// 			int vtx_ptr = k + per_proc_count_disp[cur_send_pos];
		// 			int temp_ctr = source_nodes_for_proc_count[vtx_ptr];
		// 			while(temp_ctr > 0) {
		// 				PR_accum[ outedges_in_proc[edg_ptr] - graph->start_vertex ] += PR_send[vtx_ptr];
		// 				temp_ctr--;
		// 				edg_ptr++;
		// 			}
		// 		}
		// 	}
		// 	else {
		// 		int edg_ptr = per_proc_outedges_count_disp_recv[cur_recv_pos];
		// 		for(k=0; k<per_proc_count_recv[cur_recv_pos]; k++) {
		// 			int vtx_ptr = k + per_proc_count_disp_recv[cur_recv_pos];
		// 			int temp_ctr = source_nodes_for_proc_count_recv[vtx_ptr];
		// 			while(temp_ctr > 0) {
		// 				PR_accum[ outedges_in_proc_recv[edg_ptr] - graph->start_vertex ] += PR_recv[vtx_ptr];
		// 				temp_ctr--;
		// 				edg_ptr++;
		// 			}
		// 		}
		// 	}
    //
		// 	e1 = monotonic_seconds();
		// 	int r=0;
		// 	for(r=0;r<2;r++)
		// 		MPI_Wait(&reqs[r], &status);
		// 	e2 = monotonic_seconds();
		// 	if(rank_of_the_proc==2) {
		// 		printf("Total time : %f\n", e2-s2);
		// 		printf("Computation : %f\n", e1-s1);
		// 	}
		// }

    /* Accumulate */
    double s1,s2,e1,e2,s3,e3;
    s1 = monotonic_seconds();
		// int vtx_ptr = 0;
    pr_int cnt = 0;
    pr_int* edg_ptr = outedges_in_proc_recv;
    double* vtx_ptr = PR_recv;
		for(cnt=0; cnt<total_size_recv; cnt++) {
			int temp_ctr = source_nodes_for_proc_count_recv[cnt];
      double val = *vtx_ptr;
			while(temp_ctr--) {
				PR_accum[ *edg_ptr ] += val;
        edg_ptr++;
			}
      vtx_ptr++;
		}
    // printf("Count : %d\n", cnt);
    e1 = monotonic_seconds();
    if(rank_of_the_proc == 0) {
      printf("%f\n", e1-s1);
    }


		/* Finalize new PR values */
	    double norm_changed = 0.;
	    for(v=0; v < graph->nvtxs; ++v) {
	      double const old = PR[v];
	      PR[v] = restart + (damping * PR_accum[v]);
	      norm_changed += (PR[v] - old) * (PR[v] - old);
	    }

	    double gathered_norm_change = 0.0;
	    MPI_Allreduce((void*)&norm_changed, (void*)&gathered_norm_change, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	 	if(i > 1 && sqrt(gathered_norm_change) < tol)
	 		break;
	 	free(PR_accum);
  	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank_of_the_proc==0) {
	  	double endTime = monotonic_seconds();
		printf("Global Execution Time : %f\n", endTime - startTime);
		printf("Pagerank execution : %f\n", endTime - startTime2);
	}

  	for(i=0; i<total_no_proc; i++) {
  		if(rank_of_the_proc==i) {
  			FILE * fout_final = NULL;
  			if(rank_of_the_proc==0)
  				fout_final = fopen(ofname, "w");
  			else
  				fout_final = fopen(ofname, "a");
  			for(v=0; v < graph->nvtxs; v++) {
	      		fprintf(fout_final, "%0.3e\n", PR[v]);
	    	}
  			fclose(fout_final);
  		}
  		MPI_Barrier(MPI_COMM_WORLD);
  	}

	MPI_Finalize();
	return 0;
}
