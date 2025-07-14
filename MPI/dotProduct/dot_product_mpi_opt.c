#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>


int main(int argc, char **argv){ 

    int my_PE_num, comm_size;
    int  N = 8;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size); 
    MPI_Comm_rank(MPI_COMM_WORLD,&my_PE_num);

    int A[N], B[N];
    int elements_per_proc = N/comm_size; 

    int recvbuf_A[elements_per_proc];
    int recvbuf_B[elements_per_proc]; 
    int local_sum = 0, global_sum = 0;

    if (N%comm_size != 0){
        if (my_PE_num == 0){
        printf("Error: N must be divisible by comm_size\n");
        MPI_Finalize();
        return EXIT_FAILURE;}

    } 

    if (my_PE_num == 0 ){
        for (int i=0;i<N;i++){
            A[i] = i+1;
            B[i] = 2*(i+1);
        }
    }
    clock_t start = clock();

    MPI_Scatter(A,elements_per_proc,MPI_INT,recvbuf_A,elements_per_proc,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Scatter(B,elements_per_proc,MPI_INT,recvbuf_B,elements_per_proc,MPI_INT,0,MPI_COMM_WORLD);

    for (int i=0;i<elements_per_proc;i++){
        local_sum += recvbuf_A[i]*recvbuf_B[i];
    }

    MPI_Reduce(&local_sum, &global_sum,1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    clock_t end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;

    if (my_PE_num == 0){
        printf("Final sum using MPI_Reduce: %d in %f \n", global_sum,time_taken);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}