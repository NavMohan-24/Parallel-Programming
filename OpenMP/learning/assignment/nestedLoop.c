#include <omp.h>
#include <stdio.h>
void report_num_threads(int level)
{
    #pragma omp single
    {
        printf("Level %d: number of threads in the team - %d, thread id = %d\n",
                  level, omp_get_num_threads(),omp_get_thread_num());
    }
 }
int main()
{
    omp_set_dynamic(1); // enable dynamic adjustment of number of threads.
    omp_set_nested(1);
   #pragma omp parallel num_threads(2)
    {
        report_num_threads(1);
        #pragma omp parallel num_threads(2)
        {
            report_num_threads(2);
            #pragma omp parallel num_threads(2)
            {
                report_num_threads(3);
            }
        }
    }
    return(0);
}
