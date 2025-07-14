/*
 Write a Parallel C program which should print the series of 2  and 4. Make sure both should be executed by different threads !
 */


#include <stdio.h>
#include <omp.h>
#include <math.h>

#define NUM_THREADS 2
#define N 10

int main(){
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel sections
	
	{
		#pragma omp section
	       	{
		for (int i = 0; i <N;i++){
		printf("thread= %d, series of 2 : %d\n", omp_get_thread_num(),(int)pow(2,i));
			}
		}

		#pragma omp section 
		{
			for (int i = 0; i < N;  i++){
				printf("thread= %d, series of 4 : %d\n", omp_get_thread_num(),(int)pow(4,i));

			}
		}
	}
	return 0;
}

