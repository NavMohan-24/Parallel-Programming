/*
- Initialize OpenMP with 4 threads.
- Declare an integer variable val and initialize it to a value of 1234.
- Print the initial value of val outside the OpenMP parallel region.
- Enter an OpenMP parallel region using the omp parallel directive, with the firstprivate clause applied to the variable val.
- Inside the parallel region, each thread should print the current value of val, increment it by 1, and then print the updated value.
- Print the final value of val outside the parallel region.
 */

#include <stdio.h>
#include <omp.h>


int main(){
	int val = 1234;

	printf("initial-value-of-the-variable is %d \n\n",val);
	
	omp_set_num_threads(4);

	#pragma omp parallel firstprivate(val)
	
	{
	int tid = omp_get_thread_num();
	printf("Current value of val is %d in thread number %d\n",val,tid);
	printf("Incremented value of val is %d in thread number %d\n",val+1,tid);
	}

printf("\n\nfinal value of val is %d",val);


}
