#include <stdio.h>
#include <time.h>

int main(){

    int  N = 8;
    int A[N], B[N];


    for (int i = 0; i < N; i++) {
        A[i] = i + 1;
        B[i] = (i + 1) * 2;
    }

    int RESULT = 0;
    clock_t start = clock();
    for (int i = 0; i < N; i++){

        RESULT += A[i]*B[i];

    }
    clock_t end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("found result: %d in %f seconds\n", RESULT,time_taken);
    return 0;

}