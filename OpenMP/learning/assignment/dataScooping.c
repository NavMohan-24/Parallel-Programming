/*
   1. Declare a shared integer variable (shared_variable) initialized to zero.
   2. Declare a private integer variable (private_variable) initialized to zero.
   3. Implement a parallel region using the OpenMP parallel directive.
   4. Use the shared clause to designate the shared variable within the parallel region.
   5. Use the private clause to designate the private variable within the parallel region.
   6. In the parallel region, create a loop using the OpenMP for directive to increment the shared and private variables independently.
   7. Print the values of the shared and private variables within the loop to observe their behavior with multiple threads.
   8. Print the final values of the shared and private variables after the parallel region to observe the impact of data scoping.	

*/

#include <stdio.h>
#include <omp.h>

int shared_var = 0;
int private_var = 0;
#define NUM_THREADS 4
#define SIZE 20


int main(){
	
	#pragma omp parallel shared(shared_var) private(private_var)
	{
		#pragma omp for
		for (int i = 0; i < SIZE; i++)
		{
		int tid = omp_get_thread_num();
		shared_var += 1;
		printf("Loop - 1, threadID: %d, SharedVar: %d \n", tid, shared_var);
		}

		#pragma omp for
                for (int i = 0; i < SIZE; i++)
                {
                int tid = omp_get_thread_num();
                private_var += 1;
                printf("Loop - 2, threadID: %d, PrivateVar: %d \n", tid,  private_var);
                }

	}
	
	printf("final value of Shared variable : %d \t final value of Private variable : %d\n", shared_var,private_var);
	return 0;
}

/*
NOTE:
- Using shared variable does not make the program safe from race conditions;:
*/

