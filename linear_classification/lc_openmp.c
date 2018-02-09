/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "omp.h"

// **************************************** Timer based methods BEGIN ****************************************

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

// **************************************** Timer based methods END... ****************************************

pthread_mutex_t mutexnum;

double **X;//[15170][DIMENSIONS];
double *Y;//[15170];
double *w;

int no_threads, iter, pointsPerThread;
int datapoints, dimensions;
int *noobColumnsMask;
int *w_index;
double *Xw;//[15170];
double *X2;
double error = 0.0, num = 0.0;

struct capsule {
	int startPosition;
	int dimension;
	double val;
};

void findNoobColumns() {
	int i,j;
	for(i=0;i<datapoints;i++) {
		for(j=0;j<dimensions;j++) {
			noobColumnsMask[j]|=(X[i][j]!=0);
		}
	}
}

void removeNoobColumns() {
	int i, cnt=0;
	for(i=0;i<dimensions;i++) {
	 	if(noobColumnsMask[i]==1) {
			w_index[cnt++]=i;
	 	}
	}
	dimensions=cnt;
}

int get(int idx) {
	return w_index[idx];
}

void calc_Xw() {
	int i,j;		
	#pragma omp parallel for private(i,j) num_threads(no_threads)
	for(i=0; i<datapoints; i++) {
		for(j=0; j<dimensions; j++) {
			Xw[i] = Xw[i] + X[i][get(j)] * w[j];
		}
	}	
}

void update_Xw(int j, double val) {
	int i=0;
	#pragma omp parallel for num_threads(no_threads)
	for(i=0;i<datapoints;i++) {
		Xw[i] = Xw[i] + X[i][get(j)] * (val - w[j]);
	}
}

double calcNumerator(int j) {
	int i=0;	
	double ret = 0.0, temp = 0.0;
	#pragma omp parallel for reduction(+: ret) private(temp) num_threads(no_threads)
	for(i=0; i<datapoints; i++) {
		temp = Y[i] - (Xw[i] - X[i][get(j)] * w[j]);
		ret = ret + X[i][get(j)] * temp;
	}	
	return ret;
}

void calcDenominator() {
	int i,j;
	for(j=0;j<dimensions;j++) {
		double ret = 0.0;
		for(i=0;i<datapoints;i++) {
			double temp = X[i][get(j)];
			ret = ret + temp * temp;
		}
		X2[j]=ret;
	}
}

void calc_Error() {
	double ret = 0.0, temp = 0.0;
	int i;
	#pragma omp parallel for reduction(+: ret) private(temp) num_threads(no_threads)
	for(i=0;i<datapoints;i++) {
		temp = Xw[i] - Y[i];
		ret = ret + temp * temp;
	}
	printf("Error : %lf\n", ret);
}

void printW() {
	int i;
	for(i=0;i<dimensions;i++)
		printf("%.8f , ", w[i]);
}

int main(int argc, char *argv[]) {
	if(argc!=5) {
		printf("Insufficient parameters\n");
	}
	else {
		// File I/O
		FILE *pointsFile = fopen( argv[1], "r" );
		FILE *labelFile = fopen( argv[2], "r" );
		if(pointsFile == 0)
            printf( "Could not open datapoints file\n" );
	    if(labelFile == 0)
	        printf( "Could not open labels file\n" );

		iter = atoi(argv[3]);
	    no_threads = atoi(argv[4]);

	   	fscanf(pointsFile, "%d %d", &datapoints, &dimensions);
	   	fscanf(labelFile, "%d", &datapoints);

	   	// Initialize
		int i,j;		
		X = (double **)malloc(datapoints * sizeof(double *));
		X[0] = (double *)malloc(sizeof(double) * dimensions * datapoints);
	    for(i = 0; i < datapoints; i++)
	        X[i] = (*X + dimensions * i);
	    for(i=0;i<datapoints;i++)
	    	for(j=0;j<dimensions;j++)
	    		fscanf(pointsFile, "%lf", &X[i][j]);
		Y = (double *)malloc(sizeof(double) * datapoints);
	    for(i=0;i<datapoints;i++)
	    	fscanf(labelFile, "%lf",&Y[i]);

		// **************************************** Timer starts ****************************************
		startTime = monotonic_seconds();

		pthread_t thread[no_threads];
		pthread_mutex_init(&mutexnum, NULL);
		int t, rc;	    

	    // Meaningful code
	    noobColumnsMask = (int *)calloc(dimensions, sizeof(int));
	    w_index = (int *)calloc(dimensions, sizeof(int));

	    /*
	    Not parallelizing this part since this would lead to frequent Cache invalidation of noobColumsnMask[]
	    */  
	    findNoobColumns();
	    removeNoobColumns();

	    /*
	     create and initialize 'w' with the effective no of dimensions
	     */
	    w = (double *)calloc(dimensions, sizeof(double));
	    int iterations;
		Xw = (double *)calloc(datapoints,sizeof(double));
		X2 = (double *)calloc(datapoints,sizeof(double));
				
		pointsPerThread = (datapoints+no_threads-1)/no_threads;

		/*
		Calculate Denominator terms
		Not parallelizing this part since this would lead to frequent Cache invalidation of X2[]
		*/		
		calcDenominator();

		/*
		Calculate Xw - parallelized
		*/
		calc_Xw();

	    /* 
	    Calculate Initial Error - parallelized
	    */	    
	    calc_Error(); 			

		/*
		Iterations
		*/
	    for(iterations=0; iterations<iter; iterations++) {
	    	int dim;	    		    		    
	    	for(dim=0; dim<dimensions; dim++) {		
	    		/*
	    		 Calculate Numerator
	    		 */
	    		num = calcNumerator(dim);

	    		/*
	    		 Get denominator
	    		 */
	    		double den = X2[dim];				
	    		
	    		/*
	    		 Update Xw
	    		 */
	    		update_Xw(dim, num/den); // TODO: Must parallelize
	    		
	    		/*
	    		Update the new wi
	    		*/
	    		w[dim] = num/den;
	    	}	  
	    	/*
			Update Error
	    	*/
	    	calc_Error();		    	
	    }

    	// File I/O
		fclose(pointsFile);
		fclose(labelFile);

		//printW();

		// free em resources
		free(X);
		free(Y);
		free(Xw);
	}

	// **************************************** Timer stops ****************************************
	endTime = monotonic_seconds();
	print_time(endTime-startTime);
	
	pthread_mutex_destroy(&mutexnum);
	pthread_exit(NULL);
	
	return 0;
}