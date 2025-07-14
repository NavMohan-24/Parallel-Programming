#include <stdio.h>
#include <omp.h>

#define N 4
#define CHUNK_SIZE 2
#define SIZE 10

int main(){
	float a[SIZE], b[SIZE], c[SIZE];

	for (int i =0; i<SIZE; i++){
		a[i]=(float) rand();
		b[i]=(float) rand();
	}

	omp_set_num_threads(N);

	#pragma omp parallel for schedule(static,CHUNK_SIZE)
	for (int i = 0; i < SIZE; i++)
	{
		int tid = omp_get_thread_num();
		    c[i] = a[i]+b[i];
		    printf("thread_id: %d, c[%d] = %f\n", tid,i,c[i]);	        
	}	
	return(0);
}
