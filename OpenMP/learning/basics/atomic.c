/*
Example of Atomic
*/


#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 10


int main(){
	double s = 0.0, s_local , a[N];
	for (int i=0; i<N; i++){
		a[i] = i+1;
	}

	#pragma omp parallel private(s_local)
	{
	#pragma omp for

		for (int i =0; i <N;i++)	
		{
		 s_local += a[i];
		}
		#pragma omp atomic
		s+= s_local;

		}
	printf(" value of sum : %f",s);
	return 0;
}
