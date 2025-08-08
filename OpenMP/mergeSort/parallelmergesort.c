#include <stdio.h>
#include <omp.h>

#define NTASKS = 4;

void parallelMergesort(int* A, int low, int high, int* B, int off){

    int len = high-low+1;

    if (len == 1){
        B[off] = A[low];
    }
    else{
        int* T = (int*) malloc(len*sizeof(int));

        mid = (low + high)//2;
        mid2 = mid - low + 1


        #pragma omp parallel for schedule(dynamic)
        for (int i = 0, i < NTASKS, i++){
            parallelMergesort(A, low, mid, T, 1)
            parllelMergesort(A, mid+1, higj)
        }

        parallelMerge(T, 1,mid2,mid2+1,len, B, off)


    }


}