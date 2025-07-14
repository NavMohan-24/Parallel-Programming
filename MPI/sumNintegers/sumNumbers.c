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

    
    int* counts_send = NULL; // used only by root procees
    int* displacements = NULL; // used only by root procees
    int* numArray = NULL; // used only by root procees

    int count_recv; 
    if (rank < numProcess-1){
        count_recv = numElements; //variable is defined differently for different ranks
    }
    else{
        count_recv = numElements+rem;
    }
    int* buffer_recv =  calloc(count_recv, sizeof(int));
    
    int local_sum = 0;
    int sum  = 0; // it doesnot hurt having a copy of sum in other process 
    
    if (rank==0){ 
        
        numArray = (int*) malloc(N*sizeof(int));
        counts_send = (int*) malloc(numProcess*sizeof(int));
        displacements = (int*) malloc(numProcess*sizeof(int));

        for (int i = 0; i < N; i++){
            numArray[i] = i+1;
        }

        for (int i = 0; i < numProcess; i++){
            counts_send[i] = numElements;
            displacements[i] = i*numElements; 
        }
        counts_send[numProcess-1] += rem; 
        
    }

    MPI_Scatterv(numArray, counts_send, displacements, MPI_INT, buffer_recv, count_recv, MPI_INT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < count_recv; i++){
        local_sum += buffer_recv[i];
    }

    MPI_Reduce(&local_sum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    
    free(buffer_recv);
    if (rank==0){
        free(numArray);
        free(counts_send);
        free(displacements);
        printf("sum of first %d integers = %d\n",N,sum);
    }

    MPI_Finalize();

    return 0;
}

