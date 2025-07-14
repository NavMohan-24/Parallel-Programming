#include <stdio.h>
#include <omp.h>

#define NUM_THREADS 4

int parallelPrefixSum(double arr[],int N, double res[]){
  int j = 0.0;
#pragma omp parallel for
  for (int i = 0; i < N; i++){
    double partial_sum = 0.0;
    for (int j = 0; j < i+1; j++){ //j is private since it is defined insied the for loop
      partial_sum += arr[j];
    }
    res[i] = partial_sum;
  }
  
}

int main(){
  omp_set_num_threads(NUM_THREADS);
  int size = 10;
  double A[size], B[size];

  for (int i = 0; i < size; i++){
    A[i] = i+1;
  }
  parallelPrefixSum(A, size, B);

  for (int i = 0; i < size; i++){
    printf( "B[%d] = %f\n", i, B[i]);
  }

}

/*
  - using firstprivate clause and defining ``partial_sum`` outside leading to accumulation of value
     -- i believe the private copy will hold the residue from previous iteration.
  - defining a variable inside the for loop will make it private initialised at 0.0.
  - but the approach is not scalable have **O(N^2)** time complexity  

 */
