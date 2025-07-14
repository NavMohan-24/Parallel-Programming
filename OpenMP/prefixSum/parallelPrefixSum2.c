/*

  Parallel Prefix Sum using Hillis & Steele Algorithm.

  - Calculates parallel prefix sum iteratively.
  - During ith iteration, j^th element of array whose prefix sum need to be found  will be summed with j-2^i th element.
  - The sum will become will become the j^th element of the auxillary array which is used in the next iteration.
  - All elements of the auxillary array at the indices j < 2^i is copied from the array before it.
  - The maximum number of iterations is given by log_2(n) where n is the number of elements in the array whose prefix sum has to be calculated  

**********psuedo-code**************

  for i = 0 to ceil(log_2(n))-1 do
   for j = 0 to n-1 do in parallell
      if j < 2^i then
        x^{i+1}_{j} = x^{i}_{j}

      else
        x^{i+1}_{j} = x^{i}_{j} + x^{i}_{j-2^{i}}


Note = x^{i} denotes the array used in ith iteration.
***************************************

- Hillis and steele is highly parallel algorithm with time complexity O(logn).
- The work complexity is O(logn+n) and sequential algorithm have only a work complexity of O(n).

A parallel algorithm is characterized by two metrics:

Work complexity (T1) : total number of primitive operations performed by all processors during its execution.

Depth (or Span, T_infty):  length of the longest chain of primitive operations that must be performed sequentially due to data dependencies.

 */



#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

#define NUM_THREADS 4

int parallelPrefixSum(double A[],int N){
  
	omp_set_num_threads(NUM_THREADS);

	int* temp = malloc(N * sizeof(int));
	  
	for (int d = 1;  d < N; d <<= 1){
	// d = 2^i and d<<=1 is bitwise shift {0001,0010,0100,1000 = 1,2,4,8}

       #pragma omp parallel for 
       for (int j = 0; j < N; j++){
	  if (j < d){
	    temp[j] = A[j];
	  }
	  
	  else{
	    temp[j] = A[j]+A[j-d];
	  }
	  
         }

	#pragma omp parallel for
        for (int i = 0; i < N; i++) {
            A[i] = temp[i];
          }
	  
         }
      
	free(temp);
}

int main(){
  int size = 8;	
  double A[size];

  for (int i = 0; i < size; i++){
    A[i] = i+1;
  }

  parallelPrefixSum(A, size);

  for (int i = 0; i < size; i++){
    printf( "A[%d] = %f\n", i, A[i]);
  }

  return 0;

}



