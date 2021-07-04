#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

typedef struct _array
{
	unsigned int Size;
	double* Values;
} array;

unsigned int someFunc(unsigned int s)
{
	return s + 1;
}

double getRandDouble(unsigned int i, unsigned int from, unsigned int to)
{
	unsigned int range = to - from;
	return (rand_r(&i) % range + from) + (double)(rand_r(&i) % 1000) / 1000;
}

double calcFunc_1(double x)
{
	return pow((x / exp(1.0)), 1.0 / 3.0);
}

double calcFunc_2(double x_1, double x_2)
{
	return log(fabs(tan(x_1 + x_2)));
}

array* splitArr(array Arr, unsigned int PartsCount)
{
	array* SplitedArr = malloc(sizeof(array) * PartsCount);
	unsigned int WindowSize = Arr.Size / PartsCount;
	for (int i = 0; i < PartsCount; i++)
	{
		unsigned int point = (i * WindowSize);
		SplitedArr[i].Size = (Arr.Size - point < WindowSize * 2)? Arr.Size - point : WindowSize;
		SplitedArr[i].Values = Arr.Values + point;
	}
	
	return SplitedArr;
}

void InsertionSort(array Arr)
{
	unsigned int k, i;
	for (i = 0; i < Arr.Size; i++)
	{
		k = i;
		while(k > 0 && Arr.Values[k] <= Arr.Values[k - 1])
		{
			double E  = Arr.Values[k];
			Arr.Values[k] = Arr.Values[k - 1];
			Arr.Values[k - 1] = E;
			
			k--;
		}
	}
}

int getNext(unsigned int* IndexArr, array* SortedArrs, unsigned int ArrsCount)
{
	int l = 0;
	while (l < ArrsCount && !(IndexArr[l] < SortedArrs[l].Size))
		l++;
		
	int IndNext = l;
	for (int i = l; i < ArrsCount; i++)
		if (IndexArr[i] < SortedArrs[i].Size)
		{
			if (SortedArrs[i].Values[IndexArr[i]] < SortedArrs[IndNext].Values[IndexArr[IndNext]])
				IndNext = i;
		}
			
	return IndNext;
}

void compareArr(array* SortedArrs, unsigned int ArrsCount, array ResArr)
{
	unsigned int IndexArr[ArrsCount];
	for (int i = 0; i < ArrsCount; i++)
		IndexArr[i] = 0;
	
	//double ValueNext = SortedArrs[0].Values[0];
	
	for (int j = 0; j < ResArr.Size; j++)
	{
		int IndNext = getNext(IndexArr, SortedArrs, ArrsCount);
		if (IndNext < ArrsCount)
		{
			ResArr.Values[j] = SortedArrs[IndNext].Values[IndexArr[IndNext]];
			IndexArr[IndNext]++;
		}
	}
}

#ifdef _OPENMP
	#include "omp.h"
#else
	int omp_get_num_procs() { return 1; }
	void omp_set_nested(int i) {return ; }
#endif

int main(int argc, char* argv[])
{
	const unsigned int A = 576;
	
	double StartTime = omp_get_wtime();

	unsigned long N;
	
	N = atoi(argv[1]);
	
	array Arr1;
	Arr1.Size = N;
	Arr1.Values = malloc(sizeof(double) * Arr1.Size);
	
	array Arr2;
	Arr2.Size = N / 2;
	Arr2.Values = malloc(sizeof(double) * Arr2.Size);
	
	array Arr2_copy;
	Arr2_copy.Size = N / 2;
	Arr2_copy.Values = malloc(sizeof(double) * Arr2_copy.Size);
	
	omp_set_nested(1);
	
	int Progress = 0;
	long double Sum = 0;
	
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			while (Progress != 7)
			{
				printf("%f\n", Progress / 7.0);
				sleep(1);
			}
		} // section
		
		#pragma omp section
		{
			for (unsigned long i = 0; i < 50; i++)
			{
				unsigned int s = i;
				int j;
				for (j = 0; j < N; j++)
					Arr1.Values[j] = ( rand_r(&s) % A ) + (double)( rand_r(&s) % 1000 ) / 1000;
		
				Progress = 1;
		
				for (j = 0; j < N / 2; j++)
				{
					Arr2.Values[j] = ( (rand_r(&s) % (A * 9) ) + A) + (double)( rand_r(&s) % 1000 ) / 1000;
					Arr2_copy.Values[j] = Arr2.Values[j];
				}
			
				Progress = 2;
	
				#pragma omp parallel for default(none) private(j) shared(Arr1)
				for(j = 0; j < Arr1.Size; j++)
					Arr1.Values[j] = calcFunc_1(Arr1.Values[j]);
		
				Progress = 3;
		
				#pragma omp parallel for default(none) private(j) shared(Arr2, Arr2_copy)
				for(j = Arr2.Size - 1; j > 0; j--)
					Arr2.Values[j] = calcFunc_2(Arr2_copy.Values[j], Arr2_copy.Values[j - 1]);
				Arr2.Values[0] = calcFunc_2(Arr2_copy.Values[0], 0);
		
				Progress = 4;
		
				#pragma omp parallel for default(none) private(j) shared(Arr2, Arr1)
				for (j  = 0; j < Arr2.Size; j++)
					Arr2.Values[j] = pow(Arr1.Values[j], Arr2.Values[j]);
		
				Progress = 5;
		
				int k = omp_get_num_procs();
		
				array* SplitedArr = splitArr(Arr2, k);
			
				#pragma omp parallel for default(none) private(j) shared(Arr2, SplitedArr, k)
				for (j = 0; j < k; j++)
					InsertionSort(SplitedArr[j]);
		
				compareArr(SplitedArr, k, Arr2_copy);
			
				Progress = 6;
			}
	
			int j;
			#pragma omp parallel for default(none) private(j) shared(Arr2_copy) reduction(+: Sum)
			for (j = 0; j < Arr2_copy.Size; j++)
			{
				if ( (long)(round(Arr2_copy.Values[j]) / round(Arr2_copy.Values[0])) % 2 == 0 )
					Sum += sin(Arr2_copy.Values[j]);
			}
			Progress = 7;
		
		} // section
	} // sections
	
	free(Arr1.Values);
	free(Arr2.Values);
	free(Arr2_copy.Values);		

	double FinishTime = omp_get_wtime();
	double TimeDelta_ms = 1000 * (FinishTime - StartTime);
	
	FILE* f = fopen("result.txt", "w");
	fprintf(f, "   -   |%ld|%Lf|%lf\n", N, Sum, TimeDelta_ms);
	fclose(f);

	return 0;
}
