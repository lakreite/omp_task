#include "gemm.h"
#include <omp.h>
#define THREAD_NUM    16

double bench_t_start, bench_t_end;

static
double rtclock()
{
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, NULL);
    if (stat != 0)
      printf ("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start()
{
    bench_t_start = rtclock ();
}

void bench_timer_stop()
{
    bench_t_end = rtclock ();
}

void bench_timer_print()
{
    printf ("Time in seconds = %0.6lf\n", bench_t_end - bench_t_start);
}


static
void init_array(
    int ni, 
    int nj, 
    int nk,
    //float *alpha,
    //float *beta,
    float C[ ni][nj],
    float A[ ni][nk],
    float B[ nk][nj])
    //*alpha = 1.5;
    //*beta = 1.2;
{
#pragma omp parallel shared(A,B,C) 
{
    //int i, j;
    //omp_set_num_threads( THREAD_NUM );
#pragma omp for  schedule(static)
    for (int i = 0; i < ni; i++)
      for (int j = 0; j < nj; j++)
        C[i][j] = (float) ((i*j+1) % ni) / ni;

#pragma omp for  schedule(static)
    for (int i = 0; i < ni; i++)
      for (int j = 0; j < nk; j++)
        A[i][j] = (float) (i*(j+1) % nk) / nk;

#pragma omp for  schedule(static)
    for (int i = 0; i < nk; i++)
      for (int j = 0; j < nj; j++)
        B[i][j] = (float) (i*(j+2) % nj) / nj;
}
}

static
void print_array(
    int ni, 
    int nj,
    float C[ ni][nj])
{
    int i, j;
    fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
    fprintf(stderr, "begin dump: %s", "C");
    for (i = 0; i < ni; i++)
        for (j = 0; j < nj; j++) {
            if ((i * ni + j) % 20 == 0) fprintf (stderr, "\n");
            fprintf (stderr, "%0.2f ", C[i][j]);
        }
    fprintf(stderr, "\nend   dump: %s\n", "C");
    fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}


static
void kernel_gemm(int ni, int nj, int nk,
    //float alpha,
    //float beta,
    float C[ ni][nj],
    float A[ ni][nk],
    float B[ nk][nj])
    //omp_set_num_threads( THREAD_NUM );
{
    //int i, j, k;
    //omp_set_num_threads( THREAD_NUM );
#pragma omp parallel shared(A, B, C) 
{
#pragma omp for schedule(static) 
    for (int i = 0; i < ni; i++) {
        
        for (int j = 0; j < nj; j++)
            C[i][j] *= 1.2;//beta;
    }
#pragma omp for schedule(static)
    for (int i = 0; i < ni; ++i)
    {
        for (int k = 0; k < nk; k++) {
            for (int j = 0; j < nj; j++)
                C[i][j] += A[i][k] * B[k][j] * 1.5;//alpha;
        }
    }
}
}

int main(int argc, char** argv)
{

    


    int ni = NI;
    int nj = NJ;
    int nk = NK;
    //float alpha;
    //float beta;
    float (*C)[ni][nj]; C = (float(*)[ni][nj])malloc ((ni) * (nj) * sizeof(float));
    float (*A)[ni][nk]; A = (float(*)[ni][nk])malloc ((ni) * (nk) * sizeof(float));
    float (*B)[nk][nj]; B = (float(*)[nk][nj])malloc ((nk) * (nj) * sizeof(float));

    init_array (ni, nj, nk, //&alpha, &beta,
         *C,
         *A,
         *B);

    bench_timer_start();

    kernel_gemm (
        ni, nj, nk,
        //alpha, beta,
        *C,
        *A,
        *B);

    bench_timer_stop();
    bench_timer_print();

    if (argc > 42 && ! strcmp(argv[0], "")) 
        print_array(ni, nj, *C);

    free((void*)C);
    free((void*)A);
    free((void*)B);

    return 0;
}
