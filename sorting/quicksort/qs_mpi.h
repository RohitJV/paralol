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
static inline double monotonic_seconds()
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}
/**
* @brief Output the seconds elapsed while execution.
*
* @param seconds Seconds spent on execution, excluding IO.
*/
static void print_time(double const seconds)
{
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
static void print_numbers(
    char const * const filename,
    uint32_t const * const numbers,
    uint32_t const nnumbers)
{
  FILE * fout;

  /* open file */
  if((fout = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "error opening '%s'\n", filename);
    abort();
  }

  /* write the header */
  fprintf(fout, "%d\n", nnumbers);

  /* write numbers to fout */
  for(uint32_t i = 0; i < nnumbers; ++i) {
    fprintf(fout, "%d\n", numbers[i]);
  }

  fclose(fout);
}

// **************************************** Helper methods END... ****************************************

int main(int argc, char *argv[]) {
	if(argc!=2) {
			printf("Insufficient parameters\n");
			return 0;
	}

	int npes, myrank;
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	printf("Rank : %d, Total : %d", myrank, npes);

	MPI_Finalize();
	return 0;
}s