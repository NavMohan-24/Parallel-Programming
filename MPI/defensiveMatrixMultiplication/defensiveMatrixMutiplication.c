#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <stdbool.h>


double* createMatrixA(int rows, int cols, bool TRANSPOSE ){

    double* matrixA = (double*) malloc(rows*cols*sizeof(double)); //array of pointers of each row
    if (matrixA == NULL) {
        perror("Memory allocation failed for rows");
        exit(EXIT_FAILURE); // Exit if allocation fails
    }

    // matrix with only odd numbers (change according to the matrix logic)
    if (TRANSPOSE){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                matrixA[i+(j*rows)] = 2.0*(i+j)+1.0; 
            }
        }
    }
    else{
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                matrixA[(i*cols)+j] = 2.0*(i+j)+1.0; 
            }
        }
    }
    return matrixA;
}


double* createMatrixB(int rows, int cols, bool TRANSPOSE){
    
    double* matrixB = (double*) malloc(rows*cols*sizeof(double)); //array of pointers of each row
    
    if (matrixB == NULL) {
        perror("Memory allocation failed for rows");
        exit(EXIT_FAILURE); // Exit if allocation fails
    }

    // matrix with only even numbers (change according to the matrix logic)
    if (TRANSPOSE){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                matrixB[i+(j*rows)] = 2.0*(i+j)+2.0;  
            }
        }
    }
    else{
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < cols; j++){
                matrixB[(i*cols)+j] = 2.0*(i+j)+2.0; 
            }
        }
    }
    return matrixB;
}

void computeScatterParams(int rows, int cols, int size,  int* counts, int* displs){

    int offset = 0;
    for (int i = 0; i < size; i++){
        
        int rowsPerProcs = rows / size + ((i < rows % size ? 1 : 0)); // distribution of rows across process.
        counts[i] = rowsPerProcs*cols; // num of elements need to be distributed across process - current implementation will not slice rows.
        displs[i] = offset; // offset provides starting index of elements which each process can accces.
        offset += counts[i];
    }
    return;
}

void computeGatherParams(int rows, int cols, int size, int* counts, int* displs) {
    int offset = 0;
    for (int i = 0; i < size; i++) {
        int rowsPerProc = rows / size + (i < rows % size ? 1 : 0);
        counts[i] = rowsPerProc * cols;
        displs[i] = offset;
        offset += counts[i];
    }
}

int main(int argc, char** argv){

    int size, rank;
    MPI_Init(&argc, &argv); //case-sensitive
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
   
    int row1 = atoi(argv[1]);
    int col1 = atoi(argv[2]);
    int row2 = atoi(argv[3]);
    int col2 = atoi(argv[4]);

    if (rank == 0){
        if (argc < 5){
            fprintf(stderr, "Requires Dimensions of both Matrices\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
        if (col1 != row2){
            fprintf(stderr, "These Matrices cannot be Multiplied\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        if (size > row1){
            fprintf(stderr, "Matrix cannot be efficiently distributed across process; reduce number of process \n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
       
    } // basic validation
   


    // pointers for all process 
    double* matrixA = NULL;
    double* matrixB = NULL;
    double* resMatrix = NULL;
    double* recvBuf = NULL;

    

    // scattering params
    int* counts = NULL;
    int* displs = NULL;

    if (rank == 0){
        matrixA = createMatrixA(row1,col1,false); // creating matrix to be scattered in root process
        counts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        computeScatterParams(row1, col1, size, counts, displs);     // computing scattering params  
        
        // retrieve all scattered process
        resMatrix = (double*) malloc(row1*col2*sizeof(double)); 
    }

    int rowsDistributed = row1/size + ((rank < row1 % size ? 1 : 0)); // calculate how many rows have distributed across process (unevenly with load balancing)
    int bufSize = rowsDistributed*col1; 
    recvBuf = (double*) calloc(bufSize,sizeof(double));
    double* localResMatrix = (double*) calloc(rowsDistributed*col2, sizeof(double));


    MPI_Scatterv(matrixA, counts, displs, MPI_DOUBLE, recvBuf,bufSize, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    if (rank == 0){
        free(matrixA); // since we only use recvBuf from this point
        //free(counts); // used only for scattering
        //free(displs); // used only for scattering
    }

    matrixB = createMatrixB(row2,col2,true); // creates transposed B matrix (vector) on all process'.

    // matrix multiplication
    for (int i = 0; i < rowsDistributed; i++){
        for (int j = 0; j < col2 ; j++){
            for (int k = 0; k < col1 ; k++){
                localResMatrix[i*col2+j] += recvBuf[i*col1+k]*matrixB[j*col1+k];
            }
        }
    }

    // computing gather params
    if (rank == 0) {
        computeGatherParams(row1, col2, size, counts, displs);
    }
    

    MPI_Gatherv(localResMatrix, rowsDistributed*col2, MPI_DOUBLE, resMatrix, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if (rank == 0){
        printf("Distributed Matrix Multiplication \n");
        for (int i = 0; i < row1; i++){
            for (int j = 0; j < col2; j++){
                printf("%8.1f ",resMatrix[i*col2+j]);
            }
            printf("\n");
        }


        // matrixA = createMatrixA(row1,col1,false);
        // double* finMatrix = calloc(row1 * col2, sizeof(double)); 

        // for (int i = 0; i < row1; i++){
        //     for (int j = 0; j < col2 ; j++){
        //         for (int k = 0; k < col1 ; k++){
        //             finMatrix[i*col2+j] += matrixA[i*col1+k]*matrixB[j*col1+k];
        //         }
        //     }
        // }

        // printf("Serial Matrix Multiplication \n");
        // for (int i = 0; i < row1; i++){
        //     for (int j = 0; j < col2; j++){
        //         printf("%6.1f ",finMatrix[i*col2+j]);
        //     }
        //     printf("\n");
        // }

        free(counts);
        free(displs);
        free(resMatrix);
    }

    free(matrixB);
    MPI_Finalize();

    return 0;

}

    

    

