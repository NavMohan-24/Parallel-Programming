
/*
 1. Implement a linear search algorithm using a function named linearSearch. This function should take an array, its size, and a target value as parameters.
2. Utilize OpenMP directives to parallelize the search loop within the linearSearch function. Specifically, use #pragma omp parallel for to allow multiple threads to search different parts of the array concurrently.
3. Use a critical section (#pragma omp critical) to ensure that the result variable is updated safely when the target value is found. The result variable should store the position of the target value in the array.
4. In the main function, initialize an array of integers and perform a linear search for a specified target value using the linearSearch function.
5. Print the result of the search, indicating whether the target value was found and its position if found.
 */

#include <stdio.h>
#include <omp.h>


#define NUM_THREADS 5
#define SIZE 100

int linearSearch(int arr[], int size, int val){
		
	int res;

	#pragma omp parallel for
       	for (int i = 0; i < size; i++)	
	{
		if (arr[i] == val){
			#pragma omp critical
			res = i;	
		}
	}
	return res;
}

int main(){

	omp_set_num_threads(NUM_THREADS);
	int array[SIZE];

	int result, target = 55;

	for( int i = 0; i < SIZE; i++){
		
		array[i] = i+3;
	}

	printf("result : %d\n", linearSearch(array, SIZE,target));

}

