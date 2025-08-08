#include <stdio.h>
#include <stdlib.h>

void merge(int* arr,int* aux, int low, int mid, int high){
    
    //copying array to auxillary array
    for (int k = low; k <= high; k++){
        aux[k] = arr[k];
    }

    int i = low; int j = mid+1;

        // Merge the two halves back into the original array
        for (int k = low; k <= high; k++) {
            if (i > mid) { // First half is exhausted
                arr[k] = aux[j++];
            } else if (j > high) { // Second half is exhausted
                arr[k] = aux[i++];
            } else if (aux[j] < aux[i]) { // Element from second half is smaller
                arr[k] = aux[j++];
            } else { // Element from first half is smaller or equal
                arr[k] = aux[i++];
            }
        }
}

void mergeSort(int* arr, int* aux, int l, int r){

    if (l>=r){
        return;
    }

    int mid = (l+r)/2;
     // Recursively sort the two halves
    mergeSort(arr,aux,l,mid);
    mergeSort(arr,aux,mid+1,r);
    // Merge the sorted halves
    merge(arr, aux,l,mid,r);
}

int main(){

  int N = 10;
  int* A = (int*) malloc(N * sizeof(int));
  int* aux = (int*) malloc(N * sizeof(int)); // Allocate auxiliary array once

  int values[] = {1, 4, 5, 2, 3, 6, 9, 8, 7,10};

  printf("Array A:\t");
  for (int i = 0; i < N; i++) {
      A[i] = values[i];
      printf("%d ",A[i]);
  }
  printf("\n");
  
  mergeSort(A, aux, 0, N-1);

  printf("Sorted Array A:\t");
  for (int i = 0; i < N; i++) {
     printf("%d ", A[i]);
  }
  printf("\n");

  free(A); 
  return 0;
}


