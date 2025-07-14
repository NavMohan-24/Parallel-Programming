/*
 1. Create an array of integers named array with a size of ARRAY_SIZE, and initialize it with values ranging from 0 to ARRAY_SIZE.
 2. Design a parallel region using OpenMP directives to distribute the array summation across multiple threads.
 Utilize the #pragma omp parallel and #pragma omp for directives for parallel execution.
 3. Each thread should compute its local sum by iterating over a subset of the array elements.
 4. Implement a critical section using the #pragma omp critical directive to safely update a shared sum variable with each thread's local sum.
 5. After the parallel region, calculate the average by dividing the total sum by the size of the array.
 6. Print the computed average.
 */


#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define ARRAY_SIZE 10
#define NUM_THREADS 5


int main(){
	omp_set_num_threads(NUM_THREADS);
	
	float A[ARRAY_SIZE]; // if we take as int,any arithmetic operations would consider only integer parts.
       	float s, avg = 0.0; 

	for(int i = 0; i < ARRAY_SIZE; i++)
	{
	  A[i] = i;

	 printf("A[%d] = %f\t",i, A[i]); 
	}

	#pragma omp parallel for private(s)	
	for (int i = 0; i < ARRAY_SIZE ; i++)
	{
	      s = A[i]/ARRAY_SIZE; 
	      int tid = omp_get_thread_num();

	      printf("tid : %d, s : %f\n", tid, s);
	      #pragma omp critical
	      avg += s;	

	}
	printf("avg of array : %f\n",avg);

}
