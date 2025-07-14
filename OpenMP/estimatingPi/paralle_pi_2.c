#include <stdio.h>
#include <omp.h>

#define NUM_POINTS 1000

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

   int istart = id*(NUM_POINTS/nthreads); int iend = (id+1)*(NUM_POINTS/nthreads);
   
   for(int i=istart;i<iend; i++){
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

  Code is parallelised using blocked distribution of for loop iterations.
  Each block is of size (NUM_POINTS/nthreads) and handled by seperate threads¸.

*/
