#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "fwSignal.h"
#include "fwBase.h"

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

void showArr(double* Arr, unsigned long N)
{
	printf("\n");
	int i;
	for(i = 0; i < N; i++){
		printf("%lf ", Arr[i]);
	}
	printf("\n");
}

int main(int argc, char* argv[])
{
	struct timeval StartTime, FinishTime;
	gettimeofday(&StartTime, NULL);

	unsigned long N;
	
	N = atoi(argv[1]);
	
	unsigned long NumThreads = atoi(argv[2]);
	fwSetNumThreads(NumThreads);
	
	double* Arr1 = malloc(sizeof(double) * N);
	double* Arr2 = malloc(sizeof(double) * N / 2);
	double* Arr3 = malloc(sizeof(double) * N);
	double* Arr4 = malloc(sizeof(double) * N);
	
	for(int j = 0; j<N; j++) Arr4[j] = 1.0 / 3;
	
	for (unsigned long i = 0; i < 50; i++)
	{			
		//  fill Arr1 ( A ), Arr2 (10 * A)
		unsigned int s = 0;
		int j;
		for (j = 0; j < N; j++)
			Arr1[j] = ( rand_r(&s) % A ) + (double)( rand_r(&s) % 1000 ) / 1000;
				
		for (j = 0; j < N / 2; j++)
			Arr2[j] = ( (rand_r(&s) % (A * 9) ) + A) + (double)( rand_r(&s) % 1000 ) / 1000;
		
		// calculate Arr1 / e     кубический корень 
		for(j = 0; j<N; j++)
			Arr3[j] = (Arr1[j] / exp(1.0));
		fwsPow_64f_A50(Arr3, Arr4, Arr1, N); //Arr1 = Arr3 ^ Arr4  ... Arr4 = (1.0 / 3)
		
		// ln ( tg(Arr2[i] + Arr2[i - 1] ) ) 
		for(j=N / 2; j > 0; j--)
			Arr3[j] = Arr2[j]+Arr2[j-1];
		Arr3[0] = Arr2[0];
		fwsTan_64f_A50(Arr3, Arr2, (N / 2)); //Arr2 = tg(Arr3)
		fwsAbs_64f(Arr2, Arr3, (N / 2));    //Arr3 = abs(Arr2)
		fwsLn_64f(Arr3, Arr2, (N / 2));     //Arr2 = ln(Arr3)
		
		
		// Arr2[i] = Arr1[i] ^ Arr2[i]  (  для половины Arr1)
		fwsPow_64f_A50(Arr1, Arr2, Arr3, N / 2); // Arr3 = Arr1 ^ Arr2
		
		// sort Arr3 Insertion sort  
		InsertionSort(Arr3, N / 2);
	}
	
	InsertionSort(Arr3, N / 2); // debug
	
	// calculate sum
	long double Sum = 0;
	int j;
	for (j = 0; j < N / 2; j++){
		if ( (long)(round(Arr3[j]) / round(Arr3[0])) % 2 == 0 )
			Sum += sin(Arr3[j]);
	}
		
	free(Arr1);
	free(Arr2);
	free(Arr3);
	free(Arr4);
	
	gettimeofday(&FinishTime, NULL);
	
	long TimeDelta_ms = 1000 * (FinishTime.tv_sec - StartTime.tv_sec) + (FinishTime.tv_usec - StartTime.tv_usec) / 1000;
	
	FILE* f = fopen("result.txt", "w");
	fprintf(f, "\n N=%ld NThread=%ld | sum=%Lf Milliseconds passsed: %ld\n", N, NumThreads, Sum, TimeDelta_ms);
	fclose(f);
	
	return 0;
}
