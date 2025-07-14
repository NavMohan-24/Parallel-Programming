#include <stdio.h>
#include <omp.h>
#include <stdlib.h>

#define NUM_THREADS 4


int sortArray(double* A, int N,double* sortA ){

  for (int i = 0; i < N; i++){
    sortA[i] = A[i];
  }

#pragma omp parallel for schedule(dynamic) 
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){

      if (sortA[i] < sortA[j]){
	double tmp = sortA[i];
	#pragma omp critical
	  sortA[i] = sortA[j];
	  sortA[j] = tmp ;
	}
    }
  }
  return 1;
}




int main(){

  int N = 10;
  omp_set_num_threads(NUM_THREADS);
  double A[N] ;
  for (int i = 0; i < N; i++){
    A[i] =  rand() % N;
    printf("A[%d] = %f\n", i, A[i]);
  }

  double sortA[N];

  sortArray(A, N, sortA);
  
  for (int i = 0; i < N; i++){
    printf("sortA[%d] = %f\n", i, sortA[i]);                                                                                 
  }

  return 1;
}

/*

  Problems with parallelising bubble sort algorithms:
  - even though it all critical resolves the race conditions, code still suffers from dependent writes
  - sortA[i] and sortA[j] might already be changed by another thread before your thread reaches the comparison or the update.
  - even though we are avoiding memory crashes, still it can have incorrect results or array will not be completely sorted.
  - bubbleSort cannot be parallelised efficiently.
  - to parallelise sorting algorithms, we can use Merge Sort algorithms.

 */
