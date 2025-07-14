#include <stdio.h>
#include <omp.h>

void bubble_sort(int arr[], int n) {
    int i, j;

    #pragma omp parallel for private(j) shared(arr)
    for (i = 0; i < n-1; i++) {
        for (j = 0; j < n-i-1; j++) {
            if (arr[j] > arr[j+1]){
 		int temp1 = arr[j];
		
		#pragma omp critical
		arr[j] = arr[j+1];
		arr[j+1] = temp1;		
            }
        }
    }
 
}

int main() {
   
    	int N = 10, arr[N];
	for (int i = 0 ; i < N; i++){
		arr[i] = rand() % N;
	}	
	
	bubble_sort(arr,N);

	for (int i = 0; i < N; i++){
		printf("sorted array element A[%d]: %d\n",i,arr[i]);
	}


    return 0;
}


