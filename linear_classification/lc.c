/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

float **X;//[15170][DIMENSIONS];
float *Y;//[15170];
float *w;

int no_threads, iter;
int datapoints, dimensions;
int *noobColumnsMask;
int *w_index;
float *Xw;//[15170];
float *X2;
//float* result;

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
	memset(Xw,0, sizeof(Xw));
	for(i=0;i<datapoints;i++) {
		for(j=0;j<dimensions;j++) {
			Xw[i] = Xw[i] + X[i][get(j)] * w[j];
		}
	}
}

void update_Xw(int j, float val) {
	int i;
	for(i=0;i<datapoints;i++) {
		Xw[i] = Xw[i] + X[i][get(j)] * (val - w[j]);
	}
}

float calcNumerator(int j) {
	float ret = 0.0, temp = 0.0;
	int i;
	for(i=0;i<datapoints;i++) {
		temp = Y[i] - (Xw[i] - X[i][get(j)] * w[j]);
		ret = ret + X[i][get(j)] * temp;
	}
	return ret;
}

void calcDenominator() {
	int i,j;
	for(j=0;j<dimensions;j++) {
		float ret = 0.0;
		for(i=0;i<datapoints;i++) {
			float temp = X[i][get(j)];
			ret = ret + temp * temp;
		}
		X2[j]=ret;
	}
}

void calc_Error() {
	float ret = 0.0;
	int i;
	for(i=0;i<datapoints;i++) {
		float temp = Xw[i] - Y[i];
		ret = ret + temp * temp;
	}
	printf("Error : %f\n", ret);
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
		X = (float **)malloc(datapoints * sizeof(float *));
		X[0] = (float *)malloc(sizeof(float) * dimensions * datapoints);
	    for(i = 0; i < datapoints; i++)
	        X[i] = (*X + dimensions * i);
	    for(i=0;i<datapoints;i++)
	    	for(j=0;j<dimensions;j++)
	    		fscanf(pointsFile, "%f", &X[i][j]);
		Y = (float *)malloc(sizeof(float) * datapoints);
	    for(i=0;i<datapoints;i++)
	    	fscanf(labelFile, "%f",&Y[i]);

		// **************************************** Timer starts ****************************************
		startTime = monotonic_seconds();

		pthread_t thread[no_threads];
		int t;	    

	    // Meaningful code
	    noobColumnsMask = (int *)calloc(dimensions, sizeof(int));
	    w_index = (int *)calloc(dimensions, sizeof(int));
	    memset(noobColumnsMask,0,sizeof(noobColumnsMask));
	    findNoobColumns(); // TODO: Can parallelize
	    removeNoobColumns(); // TODO: Can parallelize

	    // create and initialize 'w' with the effectine no of dimensions
	    w = (float *)calloc(dimensions, sizeof(float));

	    int iterations;
		Xw = (float *)malloc(sizeof(float) * datapoints);
		X2 = (float *)malloc(sizeof(float) * datapoints);
		calcDenominator(); // TODO: Can parallelize
	    calc_Xw(); // TODO: Can parallelize
	    calc_Error(); // TODO: Must parallelize
	    for(iterations=0; iterations<iter; iterations++) {
	    	int dim;
	    	float temp[dimensions];
	    	for(dim=0; dim<dimensions; dim++) {				
	    		float num = calcNumerator(dim); // TODO: Must parallelize
	    		float den = X2[dim];				
	    		update_Xw(dim, num/den); // TODO: Must parallelize
	    		w[dim] = num/den;
	    	}
	    	printf("\n");
	    	calc_Error(); // TODO: Must parallelize
	    }

    	// File I/O
		fclose(pointsFile);
		fclose(labelFile);

		printW();

		// free em resources
		free(X);
		free(Y);
		free(Xw);
	}

	// **************************************** Timer stops ****************************************
	endTime = monotonic_seconds();
	print_time(endTime-startTime);
	return 0;
}