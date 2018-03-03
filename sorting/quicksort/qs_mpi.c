#include <time.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

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
  fprintf(fout, "%d\n", nnumbers);

  /* write numbers to fout */
  uint32_t i;
  for(i = 0; i < nnumbers; ++i) {
    fprintf(fout, "%d\n", numbers[i]);
  }

  fclose(fout);
}

void swap(uint32_t* a, uint32_t* b) {
	uint32_t c =*a;
	*a = *b;
	*b = c;
}

void print_array(MPI_Comm comm, uint32_t* a, int no_of_elements_for_proc) {
	uint32_t size_of_the_comm, rank_of_the_proc;
	MPI_Comm_size(comm, &size_of_the_comm);
	MPI_Comm_rank(comm, &rank_of_the_proc);

	MPI_Barrier(comm);
	printf("\n..%d..\n", rank_of_the_proc);
	int i,j;
	for(i=0;i<size_of_the_comm;i++) {
		if(rank_of_the_proc==i) {
			for(j=0;j<no_of_elements_for_proc;j++) {
				printf("%d ", a[j]);
			}
			printf("\n....\n");
		}
		MPI_Barrier(comm);	
	}
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

// **************************************** Helper methods END... ****************************************

int pick_random_idx(int n, int rank_of_the_proc) {
	time_t t;
	srand((uint32_t) (time(&t)+rank_of_the_proc));
	return rand()%n;
}

uint32_t partition(MPI_Comm comm, uint32_t* a, int no_of_elements_for_proc) {	
	int i, j;
	uint32_t lesser_count=0;
	uint32_t size_of_the_comm, rank_of_the_proc;
	MPI_Comm_size(comm, &size_of_the_comm);
	MPI_Comm_rank(comm, &rank_of_the_proc);
	printf("\nInside func call with comm of size %d\n",size_of_the_comm);

	uint32_t rand_idx = pick_random_idx(no_of_elements_for_proc, rank_of_the_proc);
	printf("Element picked : %d\n", a[rand_idx]);

	/*
	All-to-all communication of random numbers, to select a pivot
	*/
	uint32_t selected_element = a[rand_idx];
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
		// if(rank_of_the_proc==0)
		// 	printf("%d : ", received_contenders[i]);
	}

	int pivot;
	if(size_of_the_comm % 2) {
		if(size_of_the_comm==1)
			pivot = received_contenders[size_of_the_comm/2];
		else
			pivot = (received_contenders[size_of_the_comm/2]+received_contenders[size_of_the_comm/2]) / 2;
	}		
	else
		pivot = received_contenders[size_of_the_comm/2];
	printf("PIVOT: %d\n", pivot);

	/*
	Perform partition
	*/
	int idx = -1;
	for(i=0;i<no_of_elements_for_proc;i++) {
		if(a[i] < pivot) {
			lesser_count++;
			idx++;
			swap(&a[i], &a[idx]);
		}
	}
	return lesser_count;
}

uint32_t calc_division_of_proc(int lesser_count, int greater_count, int total_no_proc) {
	double res = (double)(lesser_count)/(lesser_count + greater_count) * total_no_proc;
	uint32_t ret;
	if(res<1) 
		ret = _ceil(res);
	else
		ret = _floor(res);
	return ret;
}

void quicksort(MPI_Comm comm, uint32_t* a, int no_of_elements_for_proc) {	
	uint32_t size_of_the_comm, rank_of_the_proc;
	MPI_Comm_size(comm, &size_of_the_comm);
	MPI_Comm_rank(comm, &rank_of_the_proc);
	printf("\nInside func call with comm of size %d\n",size_of_the_comm);

	/*
	Calculate lower and greater than pivot elements count and broadcast it to everyone
	*/
	uint32_t lesser_count = partition(comm, a, no_of_elements_for_proc);
	uint32_t greater_count = no_of_elements_for_proc - lesser_count;
	print_array(comm, a, no_of_elements_for_proc);

	int lesser_count_cumulative, greater_count_cumulative, lesser_count_cumulative_exclusive, greater_count_cumulative_exclusive;
	MPI_Scan( (void*)&lesser_count, (void*)&lesser_count_cumulative, 1, MPI_UNSIGNED, MPI_SUM, comm);
	MPI_Scan( (void*)&greater_count, (void*)&greater_count_cumulative, 1, MPI_UNSIGNED, MPI_SUM, comm);
	lesser_count_cumulative_exclusive = lesser_count_cumulative - lesser_count;
	greater_count_cumulative_exclusive = greater_count_cumulative - greater_count;
	
	uint32_t *counts = (int *)malloc(sizeof(uint32_t) * 2);
	if(rank_of_the_proc == size_of_the_comm-1) {	
		counts[0] = lesser_count_cumulative;
		counts[1] = greater_count_cumulative;
	}
	MPI_Bcast( (void*)counts, 2, MPI_UNSIGNED, size_of_the_comm-1, comm );
	printf("Processor rank %d : both : %d, %d\n", rank_of_the_proc, counts[0], counts[1]);

	/*
	Determine the number of processors required for next iteration
	*/
	int lesser_no_proc, higher_no_proc;
	lesser_no_proc = calc_division_of_proc(counts[0], counts[1], size_of_the_comm);
	higher_no_proc = size_of_the_comm - lesser_no_proc;
	printf("Lesser proc count: %d\n", lesser_no_proc);	
	
	/*
	Determine which elements should go to which processor
	*/
}

int main(int argc, char *argv[]) {
	if(argc!=2) {
			printf("Insufficient parameters\n");
			return 0;
	}
	
	MPI_Init(&argc, &argv);

	int total_no_proc, rank_of_the_proc, total_no_of_elements, no_of_elements_per_proc;
	
	MPI_Comm_size(MPI_COMM_WORLD, &total_no_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_of_the_proc);
	printf("Rank : %d, Total : %d\n", rank_of_the_proc, total_no_proc);

	total_no_of_elements = atoi(argv[1]);
	no_of_elements_per_proc = total_no_of_elements / total_no_proc;

	// Initialize the arrays for each processor
	uint32_t *a = (uint32_t *)malloc(sizeof(uint32_t) * no_of_elements_per_proc);
	int i;
	srand((uint32_t) rank_of_the_proc+1);
	printf("Processor Rank : %d\n",rank_of_the_proc);
	for(i=0;i<no_of_elements_per_proc;i++) {
		a[i] = rand();
		printf("%d ", a[i]);
	}

	quicksort(MPI_COMM_WORLD, a, no_of_elements_per_proc);

	MPI_Finalize();
	return 0;
}