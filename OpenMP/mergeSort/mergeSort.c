#include <stdio.h>
#include <stdlib.h>



void merge(int* arr, int* L, int nL, int* R, int nR){
    
  int i = 0; int j = 0; int k = 0;

  while ((i < nL) && (j< nR)){
    if (L[i]<=R[j]){
      arr[k] = L[i];
      i += 1;
    }
    else{
      arr[k] = R[j];
      j += 1;
    }
    k += 1;
  }
  while (i < nL){
    arr[k] = L[i];
    i += 1;
    k += 1;
  }
  while (j < nR){
    arr[k] = R[j];
    j += 1;
    k += 1;
  }
}


void mergeSort(int* arr, int N){

  if (N <= 1)
    return;

  int nL = N/2;
  int nR = N - nL;

  int* L = (int*) malloc(nL*sizeof(int));
  int* R = (int*) malloc(nR*sizeof(int));

  for (int i =0 ; i < nL; i++){
    L[i] = arr[i];
  }
  for (int i = 0; i < nR; i++){
    R[i] = arr[i+nL];
  }

  mergeSort(L, nL);
  mergeSort(R, nR);
  merge(arr, L, nL, R, nR);

  free(L);
  free(R);
}


int main(){

  int N = 10;
  int* A = (int*) malloc(N * sizeof(int));

  int values[] = {1, 4, 5, 2, 3, 6, 9, 8, 7,10};

  printf("Array A:\t");
  for (int i = 0; i < N; i++) {
      A[i] = values[i];
      printf("%d ",A[i]);
  }
  printf("\n");
  
  mergeSort(A, N);

  printf("Sorted Array A:\t");
  for (int i = 0; i < N; i++) {
     printf("%d ", A[i]);
  }
  printf("\n");

  free(A); 
  return 0;
}
