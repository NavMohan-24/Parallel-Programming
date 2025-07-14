#include <stdio.h>
#include <omp.h>


#define NUM_STEPS 1000

float fun(float x){
  return 4.0/(1.0+x*x);
}

float main(){
  float h = 1.0/NUM_STEPS;
  float sum = 0.0;
  double st = omp_get_wtime();
  omp_set_num_threads(4); // if unset by default utilize complete threads.
    
#pragma omp parallel for reduction(+:sum)
  for (int i=0; i<NUM_STEPS; i++){
    
    float x_i = i*h;
    float x_j = (i+1)*h;

    //using reduction of variable     
    sum += (h/2)*(fun(x_i)+fun(x_j));

    //using critical directive
    // float res = (h/2)*(fun(x_i)+fun(x_j));
    // #pragma omp critical
    // sum += res;
  }
  double et = omp_get_wtime()-st;

  printf("result : %f\n", sum);
  printf("execution_time : %f\n", et); 
  return sum;
}
