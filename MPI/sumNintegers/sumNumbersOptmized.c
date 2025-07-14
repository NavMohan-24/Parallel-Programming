/*
MPI program to find sum of first N integers using any given number of processes. 
Example, N=10,000 and no. of processes can be 4 or 8 or 12 etc.
*/

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);
    
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // if either N or numProcess is not specified;
    if (argc < 3) {
        if (rank == 0) { // Only print error message from rank 0
            fprintf(stderr, "Usage: %s <N> <numProcess>\n", argv[0]);
        }
        MPI_Finalize();
        return 1; // Indicate error
    }
    
    int N = atoi(argv[1]);
    int numProcess = atoi(argv[2]);

    if (numProcess != size) {
        if (rank == 0) {
            fprintf(stderr, "Error: numProcess (%d) from argument does not match actual MPI size (%d).\n", numProcess, size);
        }
        MPI_Finalize();
        return 1;
    }

    
    int numElements = N/numProcess;
    int rem = N%numProcess;
    
    int local_sum = 0;
    int sum  = 0; 
    
    if (rank == numProcess-1){
        int i = rank*numElements+1;
        //printf("process %d, i value %d \n",rank,i);
        int n = i+numElements+rem;
        //printf("process %d, n value %d \n",rank,n);
        for (; i<n; i++){ //since i is already initialized
            local_sum+=i;
        }
        //printf("process %d, local_sum value %d \n",rank,local_sum);
    }
    else{
        int i = rank*numElements+1;
        //printf("process %d, i value %d \n",rank,i);
        int n = i+numElements;
        //printf("process %d, n value %d \n",rank,n);
        for (; i < n; i++){
            local_sum +=i;
        }
        //printf("process %d, local_sum value %d \n",rank,local_sum);
    }

    MPI_Reduce(&local_sum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    
    if (rank==0){
        printf("sum of first %d integers = %d\n",N,sum);
    }

    MPI_Finalize();

    return 0;
}

