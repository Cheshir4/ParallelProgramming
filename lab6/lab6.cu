#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define A 480

__global__
void map_m1(int n, float* m1) {
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= n) {
		return;
	}
	m1[i] = 1 / tanh(sqrt(m1[i]));
}

__global__
void map_m2(int n, float* m2, float* m2_copy) {
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= n) {
		return;
	}
	
	if (i != 0) {
        m2[i] += m2_copy[i];
    }

    m2[i] = fabs(tan(m2[i]));
}

__global__
void merge(int n, float* m1, float* m2) {
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= n) {
		return;
	}
	
	m2[i] = (m1[i] < m2[i]) ? m1[i] : m2[i];
}

__global__
void sort(float* m2, int start, int end) {
	int j = start;
	while (j < end - 1) {         
	    if (m2[j + 1] < m2[j]) {
	        int tmp = m2[j];
	        m2[j] = m2[j + 1];
	        m2[j + 1] = tmp;
	        j = start;
	    } else j++;
	}
}

int main(int argc, char *argv[]) {
    int N;
    float X = 0;

    struct timeval T1, T2;
    double delta_ms;
    float e = 0.00001;

    N = atoi(argv[1]); /* N равен первому параметру командной строки */
    
    gettimeofday(&T1, NULL);
    
    for(int it = 0; it < 50; it++) {
    	int j, tmp;

        unsigned int my_seed[1];
        my_seed[0] = it;
        
        /* 1. Generate */
        float *m1 = (float *) malloc(sizeof(float) * N);
        float *m2 = (float *) malloc(sizeof(float) * N / 2);
        float *m2_copy = (float *) malloc(sizeof(float) * N / 2 + 1);
        m2_copy[0] = 0;
        
        for (j = 0; j < N; j++) {
            m1[j] = rand_r(my_seed) % A + 1;
        }
        for (j = 0; j < N / 2; j++) {
            m2[j] = rand_r(my_seed) % (9 * A + 1) + A;
            m2_copy[j + 1] = m2[j];
        }
        
        /* 2. Map */
        /* Решить поставленную задачу, заполнить массив с результатами*/
        float* cm1, *cm2, *cm2_copy;
        cudaMalloc(&cm1, sizeof(float) * N);
        cudaMalloc(&cm2, sizeof(float) * N / 2);
        cudaMalloc(&cm2_copy, sizeof(float) * N / 2 + 1);

		cudaMemcpy(cm1, m1, sizeof(float) * N, cudaMemcpyHostToDevice);
		cudaMemcpy(cm2, m2, sizeof(float) * N / 2, cudaMemcpyHostToDevice);
		cudaMemcpy(cm2_copy, m2_copy, sizeof(float) * N / 2 +1, cudaMemcpyHostToDevice); 
		
		map_m1<<<(N+255)/256, 256>>>(N, cm1);  
		map_m2<<<(N+255)/256, 256>>>(N / 2, cm2, cm2_copy);  
		
        /* 3. Merge */
		merge<<<(N+255)/256, 256>>>(N / 2, cm1, cm2);  

		/* 4. Sort */
        /* Отсортировать массив с результатами указанным методом */
        cudaStream_t s;
        int step = 10;
        for(j = 0; j < N/2; j += step) {
        	cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking);
        	//printf("sorting [%d; %d]\n", j, min(N/2, j+step));
        	sort<<<1, 1, 0, s>>>(cm2, j, min(N/2, j+step));
        	cudaStreamDestroy(s);
        }
        
        cudaDeviceSynchronize();
		cudaMemcpy(m2, cm2, sizeof(float) * N / 2, cudaMemcpyDeviceToHost);
        
        j = 0;
        while (j < (N / 2) - 1) {
            if (m2[j + 1] < m2[j]) {
                tmp = m2[j];
                m2[j] = m2[j + 1];
                m2[j + 1] = tmp;
                j = 0;
            } else j++;
        }
        
        j = 0;
        while (j < (N/2) - 1)
        {
            if (m2[j+1] < m2[j])
            {
                tmp = m2[j];
                m2[j] = m2[j+1];
                m2[j+1] = tmp;
                j = 0;
            }
            else j++;
        }
        
        /* 5. Reduce */

        float min = 0;
        j = 0;
        tmp = 0;

        while (min == 0) {
            min = (fabs(m2[j]) < e) ? m2[j] : 0;
            j++;
        }
        
        for (j = 0; j < N / 2; j++) {
            int temp = 0;

            temp = (int) (m2[j] / min);

            if (temp % 2 == 0) {
                X = X + sin(m2[j]);
            }
        }
        
        cudaFree(cm1);
        cudaFree(cm2);
        cudaFree(cm2_copy);
        free(m1);
        free(m2);
        free(m2_copy);
    }
    gettimeofday(&T2, NULL);   /* запомнить текущее время T2 */
    delta_ms =  1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
    printf("\nN=%d. Milliseconds passed: %lf\n", N, delta_ms); /* T2 -T1 */
    printf("\nX: %f\n", X);
    return 0;
}
