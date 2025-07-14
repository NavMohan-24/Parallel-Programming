#include <stdio.h>
#include <omp.h>

int main() {
    int val = 42;

    #pragma omp parallel private(val)
    {
        // val is uninitialized here
        int tid = omp_get_thread_num();
        printf("Thread %d: val = %d\n", tid, val);
    }

    return 0;
}


// even though val is initialized to 42
// the val get uninitialized while entering the parallel region as private.
