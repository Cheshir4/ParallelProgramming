#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <fwBase.h>
#include <fwImage.h>
#include <fwSignal.h>

#define A 480

int main(int argc, char* argv[])
{
    int i, N, j, tmp, M;
    float X = 0;
    struct timeval T1, T2;
    long delta_ms;
    float e = 0.00001;

    N = atoi(argv[1]); /* N равен первому параметру командной строки */
    M = atoi(argv[2]); /* M равен второму параметру командной строки */
    
    //fwStaticInit();

    fwSetNumThreads(M);

    gettimeofday(&T1, NULL); /* запомнить текущее время T1 */

    for (i = 0; i < 50; i++) /* 50 экспериментов */
    {

        //srand(i);  /* инициализировать начальное значение ГСЧ   */
        unsigned int my_seed[1];
        my_seed[0] = i;

        /* 1. Generate */
        //unsigned int *my_seed = malloc(sizeof(unsigned int));
        /*float m1[N];
        float m2[N/2];
        float m2_copy[N/2+1];

        float temp[N];*/
        
        float *m1 = (float*)malloc(sizeof(float)*N);
        float *m2 = (float*)malloc(sizeof(float)*N/2);
        float *m2_copy = (float*)malloc(sizeof(float)*N/2+1);
        float *temp = (float*)malloc(sizeof(float)*N);
        
        

        for (j = 0; j < N; j++)
        {
            m1[j] = rand_r(my_seed) % A + 1;
            
            //printf("m1[%d]: %f \n", j, m1[j]);
        }

        m2_copy[0] = 0;
        for (j = 0; j < N/2; j++)
        {
            m2[j] = rand_r(my_seed) % (9*A + 1) + A;
            //printf("m2[%d]: %f \n", j, m2[j]);
            /* Подготавливаем копию для удобства сложения */
            m2_copy[j+1] = m2[j];
        }

        /* 2. Map */

        /* Работаем с первым массивом */
        /* m1[j] = 1 / tanh(sqrt(m1[j])); */
        /* Сначала вычисляем корень */
        fwsSqrt_32f_I(m1, N);
        /* Затем вычисляем гиперболический тангенс */
        fwsTanh_32f_A21(m1, temp, N);
        /* Превращаем тангенс в котангес возведением в степень -1 */
        fwsPowx_32f_A21(temp, -1, m1, N);


        /* Работаем со вторым массивом */
        /* Складываем элементы m2 и m2_copy */
        fwsAdd_32f_I(m2_copy, m2, N/2);
        /* Считаем тангенс */
        fwsTan_32f_A21(m2, temp, N/2);
        /* Считаем модуль тангенса */
        fwsAbs_32f(temp, m2, N/2);

	 //printf("blablabla\n"); /* T2 -T1 */
	
        /* 3. Merge */
        for (j=0; j<N/2; j++)
        {
            m2[j] = (m1[j] < m2[j]) ? m1[j] : m2[j];
        }

        /* 4. Sort */

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

        while (min == 0)
        {
            min = (fabs(m2[j]) < e) ? m2[j] : 0;
            j++;
        }

        for (j = 0; j < N/2; j++)
        {
            tmp = (int)(m2[j] / min);

            if (tmp % 2 == 0)
            {
                X = X + sin(m2[j]);
            }
        }
        
        free(m1);
        free(m2);
        free(m2_copy);
        free(temp);

    }
    
    gettimeofday(&T2, NULL);   /* запомнить текущее время T2 */
    delta_ms =  1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
    printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms); /* T2 -T1 */
    printf("\nX: %f\n", X); /* T2 -T1 */
    
    
    return 0;
}
