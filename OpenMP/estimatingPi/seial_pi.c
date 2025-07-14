#include <stdio.h>
#include <omp.h>

int  NUM_POINTS = 1000;



double foo(double x){
  return 4.0/(1.0+x*x);
}

int main(){

  double h = 1.0/NUM_POINTS;
  double res = 0.0;

  double st = omp_get_wtime();
  for (double i=0;i<1.0; i+=h){
    res += foo(i)*h;
    }
  double et  = omp_get_wtime();
  printf("result: %f\n",res);
  printf("time for execution: %f\n", et-st);
  return 0;
}

