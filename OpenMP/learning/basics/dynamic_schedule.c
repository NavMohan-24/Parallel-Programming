/*
 - Set the number of OpenMP threads to 2.
 - Execute a loop from 0 to 9 using OpenMP parallelization and dynamic scheduling, without specifying a chunk size.
 - Inside the loop, each thread should print its ID and the current iteration it processes.
 - Repeat the loop execution, but this time specify a chunk size of 2 for dynamic scheduling.
 - After each loop execution, print a message indicating the type of scheduling used and the results obtained. 
*/


#include <stdio.h>
#include <omp.h>

#define N 2
#define CHUNK 2

int main(){
	
	omp_set_num_threads(2);
	
	double st1 = omp_get_wtime();	
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < 10; i++){
		int tid = omp_get_thread_num();

		printf("case-1 thread_id: %d, iteration_num: %d\n", tid, i);
	}
	printf("time for completion %f",omp_get_wtime()-st1);
	
	double st2 = omp_get_wtime();
	#pragma omp parallel for schedule(dynamic,CHUNK)
	for (int i = 0; i < 10; i++){
		int tid = omp_get_thread_num();
		printf("case-2 thread_id: %d, iteration_num: %d\n", tid, i);
	}
	printf("time for completion %f", omp_get_wtime()-st2);
}





