#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>



int main(int argc, char **argv){ 
  //argument count - number of commandline arguments passed to program &
  //argument vector - array of string where each element contains one CLI argument

  int my_PE_num, comm_size;

  int  N = 8;
  MPI_Init(&argc, &argv);
  // & is the address of operator
  // passes memory address of argc and argv so that MPI could directly update
  // without &, it will pass copies of argc, argv
  
  MPI_Comm_size(MPI_COMM_WORLD,&comm_size); 
  // comm_size indicates the number of process
  // comm_size variable is populated by mpirun -n "value"
  // MPI_COMM_WORLD - denotes the default communicator of all MPI process
  
  MPI_Comm_rank(MPI_COMM_WORLD,&my_PE_num);
  int A[N], B[N], C[N]; // three arrays 
  int elements_per_proc = N/comm_size; // elements per processes
  int recvbuf_A[elements_per_proc]; // buffer to store values of A
  int recvbuf_B[elements_per_proc]; // buffer to store values of B
  int result_buf[elements_per_proc]; // buffer to store values of product

  
  if (my_PE_num==0){ // using root process to initialise A and B
    for (int i = 0; i < N; i++) {
            A[i] = i + 1;
            B[i] = (i + 1) * 2;
        }
  }

  clock_t start = clock(); 
  //starting the clock

  MPI_Scatter(A,elements_per_proc,MPI_INT,recvbuf_A,elements_per_proc,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Scatter(B,elements_per_proc,MPI_INT,recvbuf_B,elements_per_proc,MPI_INT,0,MPI_COMM_WORLD);

  // use MPI_SCATTER to distribute A(B) to recvbuf_A(B) in size of num_procs
  // num_procs denotes the number of elements are send from root process
  // MPI_INT denotes that INT data type is send
  // "0" denotes the rank of the root process.
  // MPI_COMM_WORLD denotes the MPI COMMUNICATOR

  for (int i = 0; i < elements_per_proc;i++){
    result_buf[i] = recvbuf_A[i] * recvbuf_B[i];
    // calculating product of recvbuf_A and recvbuf_B
  }

  MPI_Gather(result_buf,elements_per_proc,MPI_INT,C,elements_per_proc,MPI_INT,0,MPI_COMM_WORLD);
  //Gather data from result but to the result_buf
  
  int SUM = 0; // in the root process, the sum is calculated.
  if (my_PE_num==0){
    for (int i = 0; i < N; i++){
      SUM += C[i];
    }
    clock_t end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Final sum: %d in %f", SUM,time_taken);
  }
  MPI_Finalize();
 
  return EXIT_SUCCESS;
}


/* 

MPI_Scatter : 
- Dispatches data from the root process in the argument (rank which is specified)
- It consist of n-1 pairs of sends and receive.
- root 0 process also retains some elements for it do computation in the next step.abort
- In the example, there are 8 array and elements and if you choose number of mpi processes = 4, 
- then each process will have two elements.


*/