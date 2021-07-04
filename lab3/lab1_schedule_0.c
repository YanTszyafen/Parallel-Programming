#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>

const unsigned int A = 576;

void InsertionSort(double* Arr, unsigned long N)
{
	unsigned int k, i;
	for (i = 0; i < N; i++)
	{
		k = i;
		while(k > 0 && Arr[k] <= Arr[k - 1])
		{
			float E  = Arr[k];
			Arr[k] = Arr[k - 1];
			Arr[k - 1] = E;
			
			k--;
		}
	}
}

int main(int argc, char* argv[])
{
	struct timeval StartTime, FinishTime;
	gettimeofday(&StartTime, NULL);

	unsigned long N;
	
	N = atoi(argv[1]);
	
	double* Arr1 = malloc(sizeof(double) * N);
	double* Arr2 = malloc(sizeof(double) * N / 2);
	double* Arr2_copy = malloc(sizeof(double) * N / 2);
	
	for (unsigned long i = 0; i < 50; i++)
	{
		unsigned int s = i;
		int j;
		//  fill Arr1 ( A ), Arr2 (10 * A)
		for (j = 0; j < N; j++)
			Arr1[j] = ( rand_r(&s) % A ) + (double)( rand_r(&s) % 1000 ) / 1000;
			
		for (j = 0; j < N / 2; j++)
		{
			Arr2[j] = ( (rand_r(&s) % (A * 9) ) + A) + (double)( rand_r(&s) % 1000 ) / 1000;
			Arr2_copy[j] = Arr2[j];
		}
		
		// calculate Arr1      / e     кубический корень 
		#pragma omp parallel for default(none) private(j) shared(Arr1, N) schedule(static)
		for(j = 0; j < N; j++)
		{
			Arr1[j] = pow((Arr1[j] / exp(1.0)), 1.0 / 3.0);
		}
		
		// ln ( tg(Arr2[i] + Arr2[i - 1] ) ) 
		#pragma omp parallel for default(none) private(j) shared(Arr2, Arr2_copy, N) schedule(static)
		for(j = N / 2 - 1; j > 0; j--)
		{
			Arr2[j] = log(fabs(tan(Arr2_copy[j]+Arr2_copy[j-1])));
		}
		Arr2[0] = log(fabs(tan(Arr2[0])));
			
		// Arr2[i] = Arr1[i] ^ Arr2[i]  (  для половины Arr1)
		#pragma omp parallel for default(none) private(j) shared(Arr2, Arr1, N) schedule(static)
		for (j  = 0; j < N / 2; j++)
		{
			Arr2[j] = pow(Arr1[j], Arr2[j]);
		}
		
		// sort Arr2 Insertion sort  
		InsertionSort(Arr2, N / 2);
	}
	
	// calculate sum
	long double Sum = 0;
	int j;
	#pragma omp parallel for default(none) private(j) shared(Arr2, N, Sum) schedule(static)
	for (j = 0; j < N / 2; j++){
		if ( (long)(round(Arr2[j]) / round(Arr2[0])) % 2 == 0 )
			#pragma omp critical
			Sum += sin(Arr2[j]);
	}
	
	free(Arr1);
	free(Arr2);
	free(Arr2_copy);
	
	gettimeofday(&FinishTime, NULL);
	
	long TimeDelta_ms = 1000 * (FinishTime.tv_sec - StartTime.tv_sec) + (FinishTime.tv_usec - StartTime.tv_usec) / 1000;
	
	FILE* f = fopen("result.txt", "w");
	fprintf(f, "sheldule(static): N=%ld | sum=%Lf Milliseconds passsed: %ld\n", N, Sum, TimeDelta_ms);
	fclose(f);
	
	return 0;
}
