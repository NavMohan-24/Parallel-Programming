#include <stdio.h>
#include <omp.h>

#define NUM_POINTS 10000

double foo(double x){
  return 4.0/(1.0+x*x);}

int main(){

  double h = 1.0/NUM_POINTS;
  double res = 0.0;
  double st = omp_get_wtime();
  omp_set_num_threads(4);
  
  #pragma omp parallel reduction(+:res)
  {
   int id = omp_get_thread_num();
   int nthreads = omp_get_num_threads();
   
   for (double i=id;i<NUM_POINTS; i+=nthreads){
     double x = i*h;
     res += foo(x)*h;
    }
  }
  double et  = omp_get_wtime();
  printf("result: %f\n",res);
  printf("time for execution: %f\n", et-st);
  return 0;
}


/*

  Loop iterations are parallelised with threads.
  This type of parallelisation is called "cyclic distribution" of loop.
  Each iteration is are distributed among threads in a cyclic manner.
  Reduction clause is used to avoid **Race Conditons**.
  Race Conditions arises when different threads tries to access same resource.
  Reduction Clause creates a private copy of the **res** variable on each thread.
  Even though the variable is private copy, for reduction clause variable is initialized to neutral value ( 0 for +, 1 for *...)
  

*/
