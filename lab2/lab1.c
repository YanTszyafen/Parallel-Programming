#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

#include "fwBase.h"
#include "fwSignal.h"

/* Name = Basalaev Artem Alekseevic(h)
 * A = 8 * 5 * 10 = 400
 * Xx = 1 + ((A mod 47) mod B)
 * X1 = 1 + ((400 mod 47) mod 7) = 4
 * X2 = 1 + ((400 mod 47) mod 8) = 1
 * X3 = 1 + ((400 mod 47) mod 6) = 1
 * X4 = 1 + ((400 mod 47) mod 7) = 4
 */

void stupid(double *arr, unsigned int size) {
    unsigned int i = 1;
    while (i < size) {
        if (arr[i - 1] > arr[i]) {
            double t = arr[i - 1];
            arr[i - 1] = arr[i];
            arr[i] = t;
        } else
            i++;
    }
}


int main(int argc, char *argv[]) {
    struct timeval T1, T2;
    long time_ms;
    int N, num_threads = 1, A = 400;
    if (argc < 2) {
        printf("Usage: ./lab2-seq N [num_threads]\n");
        exit(EXIT_FAILURE);
    }
    N = atol(argv[1]); /* initialize the number N with the first command line argument */
    if (argc >= 3)
        num_threads = atol(argv[2]);
    fwSetNumThreads(num_threads);
    double *ones = malloc(sizeof(double) * N);
    fwsSet_64f((double) 1.0, ones, N); /* There is no function for dividing a constant by a 64-bit floating number*/
    gettimeofday(&T1, NULL); /* remember current time T1 */
    for (int i = 0; i < 50; ++i) {
        unsigned int int_rand = i; /* always use a different random number generator */
        /* ... here is the solution of the lab in accordance with the option */

        // Part 1. Generate
        double *M1 = malloc(sizeof(double) * N);
        for (int j = 0; j < N; ++j) {
            M1[j] = rand_r(&int_rand) * 1.0 / RAND_MAX * (A - 1) + 1.0;
        }
        double *M2 = malloc(sizeof(double) * N / 2);
        for (int j = 0; j < N / 2; ++j) {
            M2[j] = rand_r(&int_rand) * 1.0 / RAND_MAX * (9 * A) + A;
        }

        // Part 2. Map
        fwsSqrt_64f_I(M1, N);
        fwsTanh_64f_A53(M1, M1, N);
        fwsDiv_64f(M1, ones, M1, N);

        double *M2_2 = malloc(sizeof(double) * N / 2);
        M2_2[0] = sin(M2[0]);
        fwsAdd_64f(M2, (M2+1), (M2_2+1), N/2-1);
        fwsSin_64f_A53(M2_2, M2_2, N/2);
        fwsAbs_64f_I(M2_2, N/2);

        // Part 3. Merge
        fwsPow_64f_A53(M1, M2_2, M2_2, N/2);

        // Part 4. Sort
        stupid(M2_2, N / 2);

        // Part 5. Reduce

        int k = 0;
        while (M2_2[k] == 0) {
            k++;
        }
        double min = M2_2[k + 1];

        double sum = 0;
        for (int j = 0; j < N / 2; j++) {
            sum += ((int) (M2_2[j] / min) % 2) == 0 ? sin(M2_2[j]) : 0;
        }
    }
    //Check time
    gettimeofday(&T2, NULL); /* remember current time T2 */
    time_ms = 1000 * (T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
//    printf("\nN=/**/%d. \nBest time (ms): %ld\n", N, minimal_time_ms); /* elapsed time */
    printf("%d\n%ld\n", N, time_ms); /* elapsed time */

    return 0;
}