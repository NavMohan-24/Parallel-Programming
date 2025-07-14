/*
1. Initialize two arrays (vector1 and vector2) of size N with arbitrary values. You can modify the code to use your own data.

2. Use OpenMP directives to create a parallel region and distribute the loop iterations among multiple threads.

3. Implement a reduction operation to calculate the sum of corresponding elements from both vectors.

4. Ensure proper synchronization and avoid race conditions.

5. Print the final sum.
 */


#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>

#define ARRAY_SIZE 100
#define NUM_THREADS 5

int main(){
	srand(time(NULL));

	double A[ARRAY_SIZE], B[ARRAY_SIZE],sum =0.0;

	for (int i = 0; i < ARRAY_SIZE; i++){
		A[i] = (double) (rand() % 100);
		B[i] = (double) (rand() % 100);
	}
	
	#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < ARRAY_SIZE; i++){
		sum += A[i]+B[i];
	}
	printf("sum of two array: %f\n",sum);

	return 0;
}
