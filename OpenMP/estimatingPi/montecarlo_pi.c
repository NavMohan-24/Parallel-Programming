#include <stdio.h>
#include <omp.h>


#define NUM_POINTS 1000

double fun(double x){
  return 4.0/(1.0+x*x);
 }

int main(){

  double st = omp_get_wtime();
  omp_set_num_threads(4); // if unset by default utilize complete threads.
  double sum = 0.0;
  
  #pragma omp parallel reduction(+:sum){
    double x = rand() / RAND_MAX;
    sum += (1/NUM_POINTS)* fun(x);
  }

  double et = omp_get_wtime()-st;

  printf("result : %f\n", sum);
  printf("execution_time : %f\n", et); 
  return 0;
 }

