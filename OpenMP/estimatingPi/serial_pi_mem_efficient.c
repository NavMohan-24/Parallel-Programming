#include <stdio.h>
#include <omp.h>

int  NUM_POINTS = 1000;



float foo(float x){
  return 4.0/(1.0+x*x);
}

int main(){

  float h = 1.0/NUM_POINTS;
  float res = 0.0;

  float st = omp_get_wtime();
  for (float i=0;i<1.0; i+=h){
    res += foo(i)*h;
    }
  float et  = omp_get_wtime();
  
  printf("result: %f\n",res);
  printf("time for execution: %f\n", et-st);
  return 0;
}
