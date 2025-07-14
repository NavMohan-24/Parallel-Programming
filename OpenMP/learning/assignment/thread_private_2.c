/*
1. The program defines three global variables: a , b, and x. a and  x are marked as thread-private using the #pragma omp threadprivate directive.
2. In the first parallel region, each thread initializes its private copy of a , b, and x with values derived from its thread ID. Display the thread ID along with the values of a , b, and x for each thread.
3. Between the two parallel regions, simulate the master thread performing serial work. Print a message indicating that the master thread is doing serial work.
4. In the second parallel region, display the values of a , b, and x for each thread. Observe and explain any differences in the values of a and b between the first and second parallel regions1. The program defines three global variables: a , b, and x. a and  x are marked as thread-private using the #pragma omp threadprivate directive.
2. In the first parallel region, each thread initializes its private copy of a , b, and x with values derived from its thread ID. Display the thread ID along with the values of a , b, and x for each thread.
3. Between the two parallel regions, simulate the master thread performing serial work. Print a message indicating that the master thread is doing serial work.
4. In the second parallel region, display the values of a , b, and x for each thread. Observe and explain any differences in the values of a and b between the first and second parallel regions
 */



#include <stdio.h>
#include <omp.h>

#define NUM_THREADS 2

int a = 1000,b = 1000,x = 1000,tid;
#pragma omp threadprivate(a,x)


int main(){
	omp_set_num_threads(NUM_THREADS);

	#pragma omp parallel private(b,tid)
	{
	tid = omp_get_thread_num();
	a = b = x = tid;

	printf("thread ID: %d, a : %d, b : %d, c : %d\n", tid, a, b, x);
	}

	printf("Master is performing serial work here..\n");

	#pragma omp parallel private(tid)
	{
	tid = omp_get_thread_num();
	printf("thread ID: %d, a : %d, b : %d, c : %d\n", tid, a, b, x);
	}
}





// Each thread have only a private copy of b ==> thus value of b will be different across the two regions.
// since a and x is thread private, once the value is fixed in the first parallel region it propagates to second one.
