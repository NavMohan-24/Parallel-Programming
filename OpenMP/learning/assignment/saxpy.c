#include <stdio.h>
// Add the requried parallel pargram 
#include <omp.h>

void saxpy(float z[], float a, float x[], float y, int n) {
     //Declare required variables;
     int NUM_THREADS = 2;
     omp_set_num_threads(NUM_THREADS);
 	   
     // Calculate chunk size to evenly distribute the work among threads

     // Calculate the starting and ending indices for each thread

     
    // Perform SAXPY within the assigned range
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
        z[i] = a * x[i] + y;
    }
}

int main() {
    
    int n = 5;
    float z[n];
    float a = 2.0;
    float x[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    float y = 10.0;

    saxpy(z, a, x, y, n);

    // Print the result
    printf("Resulting array z: ");
    for (int i = 0; i < n; i++) {
        printf("%f ", z[i]);
    }
    printf("\n");

    return 0;
}
