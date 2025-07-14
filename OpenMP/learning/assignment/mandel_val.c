#include <stdio.h>
#include <omp.h>

#define MAXITER 200

int mandel_val(double x, double y, int maxiter) {
    int iter = 0;
    double real_part = x;
    double imag_part = y;

    while (iter < maxiter) {
        double real_temp = real_part * real_part - imag_part * imag_part + x;
        imag_part = 2 * real_part * imag_part + y;
        real_part = real_temp;

        if (real_part * real_part + imag_part * imag_part > 4.0) {
            break;
        }

        iter++;
    }

    return iter;
}

int main() {
    int m = 5 /* some value */;  // Set the value of m
    int n = 4/* some value */;  // Set the value of n
    int maxiter = MAXITER;
    int total_iters = 0;

    int depth[m][n];
    

    #pragma omp parallel for shared(depth)
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++) {
            double x = (double)i / (double)m;
            double y = (double)j / (double)n;
            depth[j][i] = mandel_val(x, y, maxiter);
           // total_iters += depth[j][i];
        }
    }

    // Rest of your code//
    for (int i=0;i<m;i++){
	    for (int j=0; j<n;j++){
		    total_iters += depth[j][i];
		}
	}
    printf("result %d\n",total_iters);
    return 0;
}
