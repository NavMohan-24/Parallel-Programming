/*
1. The program defines three global variables: a , b, and x. a and  x are marked as thread-private using the #pragma omp threadprivate directive.
2. In the first parallel region, each thread initializes its private copy of a , b, and x with values derived from its thread ID. Display the thread ID along with the values of a , b, and x for each thread.
3. Between the two parallel regions, simulate the master thread performing serial work. Print a message indicating that the master thread is doing serial work.
4. In the second parallel region, display the values of a , b, and x for each thread. Observe and explain any differences in the values of a and b between the first and second parallel regions

 */

#include <stdio.h>
#include <omp.h>

#define NUM_THREADS 2
static int  a,b,x;
#pragma omp threadprivate(x)


int main(){
	
	omp_set_num_threads(NUM_THREADS);

	#pragma omp parallel sections private(a,b)
	{
		#pragma omp section
		{
		a = omp_get_thread_num();
		b = omp_get_thread_num();
		x = omp_get_thread_num();
		printf("thread ID: %d, a: %d, b: %d, c: %d\n",omp_get_thread_num(), a, b, x);
		
		}
		
		#pragma omp section
                {
                a = omp_get_thread_num();
                b = omp_get_thread_num();
                x = omp_get_thread_num();
				
                printf("thread ID: %d, a: %d, b: %d, c: %d\n",omp_get_thread_num(), a, b, x);
                }

	}
	printf("Master thread is perfoming serial work here..\n");
}
