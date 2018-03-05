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


/**
* @brief Write an array of integers to a file.
*
* @param filename The name of the file to write to.
* @param numbers The array of numbers.
* @param nnumbers How many numbers to write.
*/
static void print_numbers(char const * const filename, uint32_t const * const numbers, uint32_t const nnumbers) {
  FILE * fout;

  /* open file */
  if((fout = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "error opening '%s'\n", filename);
    abort();
  }

  /* write the header */
  fprintf(fout, "%u\n", nnumbers);

  /* write numbers to fout */
  uint32_t i;
  for(i = 0; i < nnumbers; ++i) {
    fprintf(fout, "%u\n", numbers[i]);
  }

  fclose(fout);
}

void swap(uint32_t* a, uint32_t* b) {
	uint32_t c =*a;
	*a = *b;
	*b = c;
}

uint32_t _ceil(double x) {
	if((uint32_t)x * (double)(1.0) == x)
		return x;
	else
		return (uint32_t)x + 1;
}

uint32_t _floor(double x) {
	if((uint32_t)x * (double)(1.0) == x)
		return x;
	else
		return (uint32_t)x;
}

uint32_t _min(uint32_t x, uint32_t y) {
	if(x<y) return x;
	return y;
}

int cmpfunc (const void * a, const void * b) {
   return ( *(uint32_t*)a > *(uint32_t*)b ? 1 : -1);
}

int serial_partition(uint32_t *a, int lo, int hi) {
	int i, idx=lo-1;
	uint32_t pivot=a[hi];
	for(i=lo; i<hi; i++) {
		if(a[i] < pivot) {
			idx++;
			swap(&a[i], &a[idx]);
		}
	}
	swap(&a[hi], &a[idx+1]);
	return idx+1;
}

void serial_quicksort(uint32_t *a, int lo, int hi) {
	if(lo>=hi)
		return;
	int mid = serial_partition(a, lo,  hi);
	serial_quicksort(a, lo, mid-1);
	serial_quicksort(a, mid+1, hi);
}

// **************************************** Helper methods END... ****************************************

uint32_t pick_random_idx(uint32_t n, uint32_t rank_of_the_proc) {
	time_t t;
	srand((uint32_t) (time(&t)+rank_of_the_proc));
	return rand()%n;
}

uint32_t partition(MPI_Comm comm, uint32_t* a, uint32_t no_of_elements_for_proc) {
	uint32_t i, j, lesser_count=0, size_of_the_comm, rank_of_the_proc, pivot;
	MPI_Comm_size(comm, &size_of_the_comm);
	MPI_Comm_rank(comm, &rank_of_the_proc);

	uint32_t rand_idx = pick_random_idx(no_of_elements_for_proc, rank_of_the_proc);
	uint32_t selected_element = a[rand_idx];
	
	/*
	All-to-all communication of random numbers, to select a pivot
	*/
	uint32_t *received_contenders = malloc(sizeof(uint32_t) * size_of_the_comm);
	MPI_Allgather( (void*)&selected_element, 1, MPI_UNSIGNED, (void*)received_contenders, 1, MPI_UNSIGNED, comm);

	/*
	Sort and pick median
	*/
	for(i=0;i<size_of_the_comm;i++) {
		for(j=i+1;j<size_of_the_comm;j++) {
			if(received_contenders[i]>received_contenders[j])
				swap(&received_contenders[i],&received_contenders[j]);
		}
	}
	if(size_of_the_comm % 2)
		pivot = received_contenders[size_of_the_comm/2];
	else
		pivot = (received_contenders[size_of_the_comm/2] + received_contenders[size_of_the_comm/2 - 1]) / 2;

	/*
	Perform partition
	*/
	uint32_t idx = -1;
	for(i=0;i<no_of_elements_for_proc;i++) {
		if(a[i] < pivot) {
			lesser_count++;
			idx++;
			swap(&a[i], &a[idx]);
		}
	}
	return lesser_count;
}

uint32_t calc_division_of_proc(uint32_t lesser_count, uint32_t greater_count, uint32_t total_no_proc) {
	double res = (double)(lesser_count)/(lesser_count + greater_count) * total_no_proc;
	uint32_t ret;
	if(res<1)
		ret = _ceil(res);
	else
		ret = _floor(res);
	return ret;
}

double quicksort(MPI_Comm comm, uint32_t* a, uint32_t no_of_elements_for_proc, uint32_t total_no_proc, uint32_t total_element_count, char* filename) {
	uint32_t size_of_the_comm, rank_of_the_proc;
	MPI_Comm_size(comm, &size_of_the_comm);
	MPI_Comm_rank(comm, &rank_of_the_proc);

	/*
	If this is the only process in the communicator, 
	1. Just perform a local quicksort
	2. Distrbute all elements equally
	3. Collect everything at the last processor in MPI_COMM_WORLD
	4. Print to output file
	5. Return the time after completion - endTime
	*/
	if(size_of_the_comm == 1) {
		//qsort(a, no_of_elements_for_proc, sizeof(uint32_t), cmpfunc);
		serial_quicksort(a, 0, no_of_elements_for_proc-1);
		uint32_t proc_size, proc_rank;
		MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

		/*
		Distribute elements evenly across all processors
		*/
		uint32_t i;
		uint32_t cumulative_no_of_elements, cumulative_no_of_elements_exclusive, no_of_elements_per_proc = total_element_count / total_no_proc;
		MPI_Scan((void*)&no_of_elements_for_proc, (void*)&cumulative_no_of_elements, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
		cumulative_no_of_elements_exclusive = cumulative_no_of_elements - no_of_elements_for_proc;
		uint32_t startPosition = 0;
		uint32_t *send_count = calloc(total_no_proc, sizeof(uint32_t));
		uint32_t *receive_count = malloc(total_no_proc * sizeof(uint32_t));
		uint32_t *send_displs = calloc(total_no_proc, sizeof(uint32_t));
		uint32_t *receive_displs = calloc(total_no_proc, sizeof(uint32_t));
		for(i=0; i<total_no_proc; i++) {
			uint32_t lower_bound = i*no_of_elements_per_proc, upper_bound = (i+1)*no_of_elements_per_proc;
			uint32_t cur_position = cumulative_no_of_elements_exclusive + startPosition;
			send_displs[i] = startPosition;
			if(startPosition<no_of_elements_for_proc && (cur_position>=lower_bound && cur_position<upper_bound))
				send_count[i] = _min(upper_bound - cur_position, cumulative_no_of_elements - cur_position);
			else
				send_count[i] = 0;
			startPosition = startPosition + send_count[i];
		}
		uint32_t total_rec_cnt = 0;
		MPI_Alltoall( (void*)send_count, 1, MPI_UNSIGNED, (void*)receive_count, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
		for(i=0; i<total_no_proc; i++) {
			total_rec_cnt = total_rec_cnt + receive_count[i];
			if(i!=0)
				receive_displs[i] = receive_displs[i-1] + receive_count[i-1];
		}
		uint32_t *result_evenly_distributed = malloc(sizeof(uint32_t) * total_rec_cnt);
		MPI_Alltoallv( (void*)a, (void*)send_count, (void*)send_displs, MPI_UNSIGNED,  (void*)result_evenly_distributed, (void*)receive_count, (void*)receive_displs, MPI_UNSIGNED, MPI_COMM_WORLD);
		free(a);
		free(send_count);
		free(send_displs);
		free(receive_count);
		free(receive_displs);

		/*
		Step 2 of the algorithm completed.
		Calculate endTime.
		*/
		double endTime = monotonic_seconds();

		/*
		1. Gather all the sorted elements in the last processor of MPI_COMM_WORLD
		2. Print to file
		*/
		uint32_t *rec_cnt, *rec_displs, *result;
		if(proc_rank == proc_size-1)  {
			rec_cnt = (uint32_t*)malloc(sizeof(uint32_t) * total_no_proc);
			rec_displs = (uint32_t*)malloc(sizeof(uint32_t) * total_no_proc);		
			for(i=0; i<total_no_proc; i++) {
				rec_cnt[i] = no_of_elements_per_proc;
				rec_displs[i] = no_of_elements_per_proc*i;
			}
			result = (uint32_t*)malloc(sizeof(uint32_t) * total_element_count);
		}
		MPI_Gatherv( (void*)result_evenly_distributed, no_of_elements_per_proc, MPI_UNSIGNED, (void*)result, rec_cnt, rec_displs, MPI_UNSIGNED, proc_size-1, MPI_COMM_WORLD);

		/*
		Check correctness
		*/
		#if TEST
			if(proc_rank == proc_size-1) {
				for(i=1;i<total_element_count;i++) {
					if(result[i]<result[i-1])
						printf("NOT Sorted\n");
				}
			}
		#endif

		if(proc_rank == proc_size-1) {
			print_numbers(filename, result, total_element_count);
			free(rec_cnt);
			free(rec_displs);
			free(result);
		}
		return endTime;
	}

	/*
	1. Perform Partition
	2. Calculate lower and greater than pivot elements (cumulative) using Prefix Scan
	*/
	uint32_t lesser_count = partition(comm, a, no_of_elements_for_proc);
	uint32_t greater_count = no_of_elements_for_proc - lesser_count;
	uint32_t actual_counts[2], cumulative_counts[2], lesser_count_cumulative, greater_count_cumulative, lesser_count_cumulative_exclusive, greater_count_cumulative_exclusive;
	actual_counts[0]=lesser_count; 
	actual_counts[1]=greater_count;
	MPI_Scan( (void*)&actual_counts, (void*)&cumulative_counts, 2, MPI_UNSIGNED, MPI_SUM, comm);
	lesser_count_cumulative = cumulative_counts[0];
	greater_count_cumulative = cumulative_counts[1];
	lesser_count_cumulative_exclusive = lesser_count_cumulative - lesser_count;
	greater_count_cumulative_exclusive = greater_count_cumulative - greater_count;

	/*
	Broadcast info on total number of lower and higher than pivot elements
	*/
	uint32_t *counts = (uint32_t *)malloc(sizeof(uint32_t) * 2);
	if(rank_of_the_proc == size_of_the_comm-1) {
		counts[0] = lesser_count_cumulative;
		counts[1] = greater_count_cumulative;
		/* If no progress, repeat again */
		if(counts[0]==0 || counts[1]==0)
			return quicksort(comm, a, no_of_elements_for_proc, total_no_proc, total_element_count, filename);
	}
	MPI_Bcast( (void*)counts, 2, MPI_UNSIGNED, size_of_the_comm-1, comm );

	/*
	Determine the number of processors required for next iteration
	*/
	uint32_t lesser_no_proc, greater_no_proc;
	lesser_no_proc = calc_division_of_proc(counts[0], counts[1], size_of_the_comm);
	greater_no_proc = size_of_the_comm - lesser_no_proc;

	/*
	Determine which elements should go to which processor
	*/
	uint32_t i,j=0;
	uint32_t *send_count = calloc(size_of_the_comm, sizeof(uint32_t));
	uint32_t *receive_count = malloc(size_of_the_comm * sizeof(uint32_t));
	uint32_t *send_displs = calloc(size_of_the_comm, sizeof(uint32_t));
	uint32_t *receive_displs = calloc(size_of_the_comm, sizeof(uint32_t));
	uint32_t lesser_per_proc_count = counts[0] / lesser_no_proc;
	for(i=0; i<lesser_count; i++) {
		uint32_t dest = (lesser_count_cumulative_exclusive + i);
		while(j!=lesser_no_proc-1 && dest >= (j+1)*lesser_per_proc_count)
			j++;
		send_count[j]++;
	}
	j=0;
	uint32_t greater_per_proc_count = counts[1] / greater_no_proc;
	for(i=0;i<greater_count;i++) {
		uint32_t dest = (greater_count_cumulative_exclusive + i);
		while(j!=greater_no_proc-1 && dest >= (j+1)*greater_per_proc_count)
			j++;
		send_count[lesser_no_proc+j]++;
	}

	/*
	Calculate the send displacement
	*/
	uint32_t snd_cnt = 0;
	for(i=0;i<size_of_the_comm;i++) {
		snd_cnt = snd_cnt + send_count[i];
		if(i!=0)
			send_displs[i] = send_displs[i-1] + send_count[i-1];
	}

	/*
	1. Use AllToAll Gather to share sendcount on receiving ends, to determine receive counts
	2. Calculate receive displacement
	*/
	uint32_t rec_cnt = 0;
	MPI_Alltoall( (void*)send_count, 1, MPI_UNSIGNED, (void*)receive_count, 1, MPI_UNSIGNED, comm);
	for(i=0;i<size_of_the_comm;i++) {
		rec_cnt = rec_cnt + receive_count[i];
		if(i!=0)
			receive_displs[i] = receive_displs[i-1] + receive_count[i-1];
	}

	/*
	Use the above info to set displs in AllToAllGatherV
	*/
	uint32_t *temp = malloc(sizeof(uint32_t) * rec_cnt);
	MPI_Alltoallv( (void*)a, (void*)send_count, (void*)send_displs, MPI_UNSIGNED,       (void*)temp, (void*)receive_count, (void*)receive_displs, MPI_UNSIGNED,    comm );
	free(a);
	free(send_count);
	free(send_displs);
	free(receive_count);
	free(receive_displs);

	/*
	Split the communicators into 2 parts
	*/
	MPI_Comm newComm;
	MPI_Comm_split(comm, (rank_of_the_proc/lesser_no_proc?1:0), 1, &newComm);
	MPI_Barrier(newComm);
	return quicksort(newComm, temp, rec_cnt, total_no_proc, total_element_count, filename);
}

int main(int argc, char *argv[]) {
	if(argc!=3) {
		printf("Incorrect number of parameters passed\n");
		return 0;
	}

	MPI_Init(&argc, &argv);
	uint32_t total_no_proc, rank_of_the_proc, total_no_of_elements, no_of_elements_per_proc;
	MPI_Comm_size(MPI_COMM_WORLD, &total_no_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_of_the_proc);

	/*
	Check the hosts on which processes run
	*/
    // char proc_name[total_no_proc];
    // uint32_t name_len;
    // MPI_Get_processor_name(proc_name, &name_len);
    // printf("Hello World from processor %d out of %d, executing on %s\n",
    //       rank_of_the_proc, total_no_proc, proc_name);

	total_no_of_elements = atoi(argv[1]);
  	char* filename = argv[2];
	no_of_elements_per_proc = total_no_of_elements / total_no_proc;

	/*
	Initialize the arrays for each processor
	*/
	uint32_t *a = (uint32_t *)malloc(sizeof(uint32_t) * no_of_elements_per_proc);
	uint32_t i;
	srand((uint32_t) rank_of_the_proc);
	for(i=0;i<no_of_elements_per_proc;i++)
		a[i] = (uint32_t)rand();

	/* Test with serial quicksort : output stored in test.txt */
	#if TEST
		MPI_Barrier(MPI_COMM_WORLD);
		uint32_t *allElements;
		if(rank_of_the_proc == total_no_proc-1)
			allElements = (uint32_t*)malloc(total_no_of_elements * sizeof(uint32_t));
		MPI_Gather((void*)a, no_of_elements_per_proc, MPI_UNSIGNED, (void*)allElements, no_of_elements_per_proc, MPI_UNSIGNED, total_no_proc-1, MPI_COMM_WORLD);
		if(rank_of_the_proc == total_no_proc-1) {
			double serialStartTime, serialEndTime;
			serialStartTime = monotonic_seconds();
			qsort(allElements, total_no_of_elements, sizeof(uint32_t), cmpfunc);
			serialEndTime = monotonic_seconds();
			printf("Serial "); 
			print_time(serialEndTime - serialStartTime);
			print_numbers("test.txt", allElements, total_no_of_elements);
			free(allElements);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	#endif

	/*
	Set start time by picking the earliest starttime of all processes
	Set end time by picking the latest endtime of all processes
	*/
	double startTime, endTime;
	startTime = monotonic_seconds();
	MPI_Scan( (void*)&startTime, (void*)&startTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	endTime = quicksort(MPI_COMM_WORLD, a, no_of_elements_per_proc, total_no_proc, total_no_of_elements, filename);
	MPI_Scan( (void*)&endTime, (void*)&endTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	if(rank_of_the_proc == total_no_proc-1)
		print_time(endTime-startTime);

	MPI_Finalize();
	return 0;
}