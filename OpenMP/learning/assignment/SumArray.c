/*
 1. Initialize an array of N numbers.
 
 2. Use OpenMP directives to create a parallel region and distribute the loop iterations among multiple threads.
 
 3. Implement a reduction operation to calculate the sum of the array elements in parallel.
 
 4. Print the final sum.
 */

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>

#define ARRAY_SIZE 100
#define NUM_THREADS 5

int main(){

        srand(time(NULL));

        double A[ARRAY_SIZE], sum = 0.0;

        for (int i = 0; i < ARRAY_SIZE; i++){
                A[i] = (double) (rand() % 100);
        }

        #pragma omp parallel for reduction(+:sum)
        for (int i = 0; i < ARRAY_SIZE; i++){
                sum += A[i];
        }

        printf("sum of two array: %f\n",sum);
        return 0;
}
