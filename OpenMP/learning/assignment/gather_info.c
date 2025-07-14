/*
1. Use OpenMP routines or functions to gather information about the OpenMP execution environment.
 - The number of processors available. 
 - The number of threads currently in use.
 - The maximum number of threads that can be utilized.
 - Determine if the execution is currently within a parallel region.
 - Check and display if dynamic threads are enabled.
 - Verify and print whether nested parallelism is supported.

2. Implement the necessary queries inside the master thread.
3. Display the gathered information clearly, indicating the purpose of each output.
4. Account for potential variations in the implementation of Fortran functions (e.g., some compilers may use integers instead of logical values)
 */


#include <stdio.h>
#include <omp.h>

#define NUM_THREADS 2

int main(){
	
	omp_set_num_threads(NUM_THREADS);
	
	#pragma omp parallel 
	{
	
		#pragma omp master
		{
		 printf("num_of_processors : %d\n", omp_get_num_procs());
		 printf("num of threads used : %d\n", omp_get_num_threads());
		 printf("max threads that can be utilized : %d\n", omp_get_max_threads());
		 printf("is it in parallel region : %d\n", omp_in_parallel());
		 printf("has dynamic threads are enabled : %d\n",  omp_get_dynamic());
		 printf("whether nested parallelism is supported: %d\n", omp_get_nested());
	
		}
	}

}
