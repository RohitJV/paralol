#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

void swap(uint32_t* a, uint32_t* b) {
    uint32_t c =*a;
    *a = *b;
    *b = c;
}

uint32_t partition(uint32_t array[], uint32_t start, uint32_t end)
{
    uint32_t pivot = array[end];
    uint32_t pIndex = start;
    
    uint32_t i;
    for(i = start; i < end; i++) {
        if(array[i] <= pivot) {
            if(pIndex != i) {
                swap(&array[i], &array[pIndex]);
            }
            pIndex++;
        }
    }
    if(pIndex != end) {
        swap(&array[pIndex], &array[end]);
    }
    return pIndex;
}

void quickSort(uint32_t array[], uint32_t start, uint32_t end)
{
    if(start < end) {
        uint32_t pIndex = partition(array, start, end);
        quickSort(array, start, pIndex-1);
        quickSort(array, pIndex+1, end);
    }
}

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




int main(int argc, char *argv[])
{
    uint32_t total_no_of_elements = atoi(argv[1]);
    uint32_t size = total_no_of_elements;
    uint32_t* array = (uint32_t*)malloc(sizeof(uint32_t) * total_no_of_elements);
    time_t t;
    srand((uint32_t) (time(&t)));
    uint32_t i;
    for(i=0;i<size;i++) {
        array[i] = rand();
    }

    printf("Here");


    double startTime, endTime;
    startTime = monotonic_seconds();

    /*
     first and last indices are passed
     idea is to move lower elements to the left of the list/pivot
     */
    quickSort(array, 0, size-1);

    endTime = monotonic_seconds();

    print_time(endTime-startTime);
    // for(i = 0; i < size; i++) {
    //     printf("%d ", array[i]);
    // }
}