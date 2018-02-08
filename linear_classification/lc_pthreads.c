/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

void* calc_Xw(void *t) {
	long startPosition = (long)t;
	int endPosition = datapoints;
   	if(startPosition+pointsPerThread < endPosition)
   		endPosition=startPosition+pointsPerThread;
	int i,j;	
	for(i=startPosition;i<endPosition;i++) {
		double temp = 0.0;
		for(j=0;j<dimensions;j++) {
			temp = temp + X[i][get(j)] * w[j];
		}
		Xw[i]=temp;
	}
	pthread_exit((void *) 0);
}

void* update_Xw(void *c) {
	int i;
   	struct capsule *tempC = (struct capsule*)c;
   	int startPosition = tempC->startPosition, j = tempC->dimension;
   	double result = 0.0, temp = 0.0, val = tempC->val;
   	//printf("Thread %ld starting...\n",startPosition);
   	int endPosition = datapoints;
   	if(startPosition+pointsPerThread < endPosition)
   		endPosition=startPosition+pointsPerThread;
   	for (i=startPosition; i<endPosition; i++) {
   		Xw[i] = Xw[i] + X[i][get(j)] * (val - w[j]);
   	}
   	free(c);
   	pthread_exit((void *) 0);
}

void* calcNumerator(void *c) {
	int i;
   	struct capsule *tempC = (struct capsule*)c;
   	int startPosition = tempC->startPosition, j = tempC->dimension;
   	double result = 0.0, temp = 0.0;   	
   	int endPosition = datapoints;
   	if(startPosition+pointsPerThread < endPosition)
   		endPosition=startPosition+pointsPerThread;
   	for (i=startPosition; i<endPosition; i++) {
   		temp = Y[i] - (Xw[i] - X[i][get(j)] * w[j]);
		result = result + X[i][get(j)] * temp;			
   	}
   	pthread_mutex_lock (&mutexnum);
   	num += result;
   	pthread_mutex_unlock (&mutexnum);
   	free(c);
   	pthread_exit((void *) 0);
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

void* calc_Error(void *t) {
	double ret = 0.0;
	long startPosition = (long)t;
   	int endPosition = datapoints, i;
   	if(startPosition+pointsPerThread < endPosition)
   		endPosition=startPosition+pointsPerThread;	
	for(i=startPosition;i<endPosition;i++) {
		double temp = Xw[i] - Y[i];
		ret = ret + temp * temp;
	}
	pthread_mutex_lock (&mutexnum); 	
	error+=ret;
	pthread_mutex_unlock (&mutexnum);
	pthread_exit((void *) 0);	
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
		for(t=0; t<no_threads; t++) {
			rc = pthread_create(&thread[t], NULL, calc_Xw, (void *)(long)(t*pointsPerThread));
			if (rc) {
	        	printf("ERROR; return code from pthread_create() is %d\n", rc);
	         	exit(-1);
	        }					        
		}
		for(t=0; t<no_threads; t++) {
			pthread_join(thread[t], NULL);
		}	    

	    /* 
	    Calculate Initial Error - parallelized
	    */	    
	    error = 0.0;
	    for(t=0; t<no_threads; t++) {
			rc = pthread_create(&thread[t], NULL, calc_Error, (void *)(long)(t*pointsPerThread));
			if (rc) {
	        	printf("ERROR; return code from pthread_create() is %d\n", rc);
	         	exit(-1);
	        }					        
		}
		for(t=0; t<no_threads; t++) {
			pthread_join(thread[t], NULL);
		}	    
		printf("Initial Error : %lf\n", error);	  			

		/*
		Iterations
		*/
	    for(iterations=0; iterations<iter; iterations++) {
	    	int dim;	    		    	
	    	for(dim=0; dim<dimensions; dim++) {		
	    		/*
	    		 Calculate Numerator
	    		 */
	    		num = 0.0;		
	    		for(t=0; t<no_threads; t++) {
	    			struct capsule *c = malloc(sizeof(struct capsule));
	    			c->dimension = dim;
	    			c->startPosition = t * pointsPerThread;
	    			rc = pthread_create(&thread[t], NULL, calcNumerator, (void *)c);
	    			if (rc) {
			        	printf("ERROR; return code from pthread_create() is %d\n", rc);
			         	exit(-1);
			        }					        
	    		}
	    		for(t=0; t<no_threads; t++) {
	    			pthread_join(thread[t], NULL);
	    		}

	    		/*
	    		 Get denominator
	    		 */
	    		double den = X2[dim];				
	    		
	    		/*
	    		 Update Xw
	    		 */
	    		//update_Xw(dim, num/den); // TODO: Must parallelize
	    		for(t=0; t<no_threads; t++) {
	    			struct capsule *c = malloc(sizeof(struct capsule));
	    			c->dimension = dim;
	    			c->startPosition = t * pointsPerThread;
	    			c->val = num/den;
	    			rc = pthread_create(&thread[t], NULL, update_Xw, (void *)c);
	    			if (rc) {
			        	printf("ERROR; return code from pthread_create() is %d\n", rc);
			         	exit(-1);
			        }					        
	    		}
	    		for(t=0; t<no_threads; t++) {
	    			pthread_join(thread[t], NULL);
	    		}
	    		
	    		/*
	    		Update the new wi
	    		*/
	    		w[dim] = num/den;
	    	}	  
	    	/*
			Update Error
	    	*/
	    	error = 0.0;
		    for(t=0; t<no_threads; t++) {
				rc = pthread_create(&thread[t], NULL, calc_Error, (void *)(long)(t*pointsPerThread));
				if (rc) {
		        	printf("ERROR; return code from pthread_create() is %d\n", rc);
		         	exit(-1);
		        }					        
			}
			for(t=0; t<no_threads; t++) {
				pthread_join(thread[t], NULL);
			}	    
			printf("Error after iteration %d : %lf\n", iterations+1, error);			
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