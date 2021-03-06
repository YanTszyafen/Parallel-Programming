#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

const unsigned int A = 576;

void InsertionSort(double* Arr, unsigned long N)
{
	unsigned int k, i;
	for (i = 0; i < N; i++){
		k = i;
		while(k > 0 && Arr[k] >= Arr[k - 1]){
			double E  = Arr[k];
			Arr[k] = Arr[k - 1];
			Arr[k - 1] = E;
			
			k--;
		}
	}
}

int main(int argc, char* argv[])
{
	unsigned long i = 0, N;
	struct timeval StartTime, FinishTime;
	long  TimeDelta_ms = 0;
	
	N = atoi(argv[1]);
	
	gettimeofday(&StartTime, NULL);
	
	double* Arr1 = malloc(sizeof(double) * N);
	
	double* Arr2 = malloc(sizeof(double) * N / 2);

	for (i = 0; i < 1; i++)
	{
		
			srand(i);
			
			//  fill Arr1 ( A ), Arr2 (10 * A)
			unsigned int s = i;
			int j;
			for (j = 0; j < N; j++)
			{
				Arr1[j] = ( rand_r(&s) % A ) + (double)( rand_r(&s) % 1000 ) / 1000;
			}
			
			for (j = 0; j < N / 2; j++)
			{
				Arr2[j] = ( (rand_r(&s) % (A * 9) ) + A) + (double)( rand_r(&s) % 1000 ) / 1000;
			}
			
			// calculate Arr1      / e     кубический корень 
			for(j = 0; j<N; j++)
			{
				Arr1[j] = pow((Arr1[j] / exp(1.0)), 1/3);
			} //*/
			
			// ln ( tg(Arr2[i] + Arr2[i - 1] ) ) 
			for(j=1; j<N/2;j++)
			{
				Arr2[j] = log(fabs(tan(Arr2[j]+Arr2[j-1])));
			} 
			
			// Arr2[i] = Arr1[i] ^ Arr2[i]  (  для половины Arr1)
			for (j  = 0; j < N / 2; j++)
			{
				Arr2[i] = pow(Arr1[j], Arr2[j]);
			}//*/
		
		// sort Arr2 Insertion sort  
		
		InsertionSort(Arr2, N / 2);
		

	}
    // calculate sum
    long double Sum = 0;
    for (j = 0; j < N / 2; j++){
        if ( (long)(round(Arr2[j]) / round(Arr2[N / 2 - 1])) % 2 == 0 )
            Sum += sin(Arr2[j] ) ;
    }//*/
    //printf("\n%Lf", Sum);
	free(Arr1);
	free(Arr2);
	
	gettimeofday(&FinishTime, NULL);
	
	TimeDelta_ms = 1000 * (FinishTime.tv_sec - StartTime.tv_sec) + (FinishTime.tv_usec - StartTime.tv_usec) / 1000;
	
	printf("\n N=%ld. Milliseconds passsed: %ld\n", N, TimeDelta_ms);
	
	return 0;
}