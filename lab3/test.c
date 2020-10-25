#include <omp.h>
#include <stdio.h>

void main()
{
	int i = 1;
	#pragma omp parallel
	{
		
		
		
		//#pragma omp atomic
		i++;
		printf("%d\n", i);
		#pragma omp parallel
		{
		#pragma omp parallel
		printf("Hello, world!\n");
		}
	}
}

