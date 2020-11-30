#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#ifdef _OPENMP
#include "omp.h"
#else

int omp_set_num_threads(int M) { return 1; }

int omp_set_nested(int l) { return 1; }

#endif

#define A 480


int main(int argc, char *argv[]) {
    int i, N, M, pr;
    float X = 0;

    //struct timeval T1, T2;
    double delta_ms;
    float e = 0.00001;

    N = atoi(argv[1]); /* N равен первому параметру командной строки */
    M = atoi(argv[2]);

#ifdef _OPENMP
    double t1, t2, time1, time2;
    t1 = omp_get_wtime();
    time1 = t1;
#else
    struct timeval T1, T2, Time1, Time2;
    gettimeofday(&T1, NULL); /* запомнить текущее время T1 */
    Time1 = T1;
#endif


    //gettimeofday(&T1, NULL); /* запомнить текущее время T1 */


    omp_set_num_threads(M);
    omp_set_nested(1);

    //printf("\nReadyyy.... GO!"); 


    pr = 0;
#pragma omp parallel shared(pr)
    {
#pragma omp sections nowait
        {

#pragma omp section
            {
                //#pragma omp parallel for default(none) private(i) shared(N,e) reduction(+:X)
                for (i = 0; i < 10; i++) /* 50 экспериментов */
                {
                
#ifndef _OPENMP
			gettimeofday(&Time2, NULL);
			delta_ms = 1000 * (Time2.tv_sec - Time1.tv_sec) + (Time2.tv_usec - Time1.tv_usec) / 1000;
                if (delta_ms >= 1000) {
                	printf("\nProgram is ready on %i %%\n", pr * 2);
                    gettimeofday(&Time1, NULL);
                }
                
             //printf("Iteration  is started\n");   

#endif 


                    //printf("Iteration %i is started\n", i);

                    int j, tmp;

                    unsigned int my_seed[1];
                    my_seed[0] = i;

                    //printf("blablabla \n");
                    /* 1. Generate */

                    float *m1 = (float *) malloc(sizeof(float) * N);
                    float *m2 = (float *) malloc(sizeof(float) * N / 2);
                    float *m2_copy = (float *) malloc(sizeof(float) * N / 2 + 1);
                    m2_copy[0] = 0;
                    //float *temp = (float*)malloc(sizeof(float)*N);

                    // printf("m1[%d]: %f \n", 0, m1[0]);
                    //#pragma omp parallel sections default(none) private(j) shared(m1, N, m2, m2_copy, my_seed)
                    //{
                    //НЕ ПАРАЛЛЕЛИТСЯ, ТАК КАК RAND_R МЕНЯЕТ ПЕРЕМЕННУЮ MY_SEED, И ПРИ ЛЮБОЙ ПОПЫТКЕ РАСПАРАЛЛЕЛИТЬ ПОЛУЧИТСЯ ОТЛИЧНЫЕ ОТ ПОСЛЕДОВАТЕЛЬНОГО ЗАПОЛНЕНИЯ РЕЗУЛЬТАТЫ В МАССИВЕ И, КАК СЛЕДСТВИЕ, НЕПРАВИЛЬНЫЙ РЕЗУЛЬТАТ.
                    //#pragma omp section
                    //{
                    for (j = 0; j < N; j++) {

                        m1[j] = rand_r(my_seed) % A + 1;

                        //printf("m1[%d]: %f \n", j, m1[j]);
                    }
                    //}


                    //#pragma omp section
                    ///{
                    //#pragma omp parallel for default(none) private(...) shared(...)
                    //  НЕ ПАРАЛЛЕЛИТСЯ ПО ТОЙ ЖЕ ПРИЧИНЕ
                    for (j = 0; j < N / 2; j++) {
                        m2[j] = rand_r(my_seed) % (9 * A + 1) + A;
                        //printf("m2[%d]: %f \n", j, m2[j]);
                        /* Подготавливаем копию для удобства сложения */
                        m2_copy[j + 1] = m2[j];
                    }
                    //}
//}
#pragma omp parallel for default(none) private(j) shared(N, m2, m2_copy)
                    for (j = 0; j < N / 2; j++) {
                        /* Подготавливаем копию для удобства сложения */
                        m2_copy[j + 1] = m2[j];
                    }





                    /* 2. Map */
                    /* Решить поставленную задачу, заполнить массив с результатами*/
#pragma omp parallel for default(none) private(j) shared(N, m1)
                    for (j = 0; j < N; j++) {
                        m1[j] = 1 / tanh(sqrt(m1[j]));
                    }

#pragma omp parallel for default(none) private(j) shared(N, m2, m2_copy)
                    for (j = 0; j < N / 2; j++) {
                        if (j != 0) {
                            m2[j] = m2[j] + m2_copy[j];
                        }

                        m2[j] = fabs(tan(m2[j]));
                    }


                    /* Этап Merge */
#pragma omp parallel for default(none) private(j) shared(N, m1, m2)
                    for (j = 0; j < N / 2; j++) {
                        m2[j] = (m1[j] < m2[j]) ? m1[j] : m2[j];
                    }

                    /* 4. Sort */

                    /* Отсортировать массив с результатами указанным методом */

                    /* j = 0;
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
                     }*/

                    /* for (int i = 0; i < N/2; i++)
                    {
                        printf("\nm[%i]: %f\n", i, m2[i]);
                    } */

                    int step = (N / 2) / M;
  
                    int test = N / 2;
                    if ((test % M) > 0) { step++; }
                    if (step == 0) { step++; }                
#pragma omp parallel for default(none) private(j, tmp, i) shared(N, m2, step)
                    for (int i = 0; i <= (N/2 - step); i += step) {
                        j = i;
                      //printf("%i\n", i); 
                        while (j < (i + step - 1)) {
                         
                            if (m2[j + 1] < m2[j]) {
                                tmp = m2[j];
                                m2[j] = m2[j + 1];
                                m2[j + 1] = tmp;
                                j = i;
                            } else j++;
                        }
                    }


                    j = 0;
                    while (j < (N / 2) - 1) {
                        if (m2[j + 1] < m2[j]) {
                            tmp = m2[j];
                            m2[j] = m2[j + 1];
                            m2[j + 1] = tmp;
                            j = 0;
                        } else j++;
                    }

                    /* Этап Reduce */

                    float min = 0;
                    j = 0;
                    tmp = 0;

                    while (min == 0) {
                        min = (fabs(m2[j]) < e) ? m2[j] : 0;
                        j++;
                    }

#pragma omp parallel for default(none) private(j) shared(N, min, m2)  reduction(+:X)
                    for (j = 0; j < N / 2; j++) {
                        int temp = 0;

                        temp = (int) (m2[j] / min);

                        if (temp % 2 == 0) {
                            X = X + sin(m2[j]);
                        }
                    }
                    pr++;
                    //printf("Iteration %i is ended\n", i);
                }
            }
            
            
#pragma omp section
            {

#ifdef _OPENMP
                while (pr < 50) {
                    time2 = omp_get_wtime();
                    if (time2 - time1 >= 1) {
                        printf("\nProgram is ready on %i %%\n", pr * 2);
                        time1 = omp_get_wtime();
                    }
                }
                 //printf("Iteration  is dfttghd\n");  
       
#endif
			}
        }
        
       
    }


#ifdef _OPENMP
    t2 = omp_get_wtime();
    delta_ms = (t2-t1)*1000;
#else
    gettimeofday(&T2, NULL);   /* запомнить текущее время T2 */
    delta_ms = 1000 * (T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
#endif
    //gettimeofday(&T2, NULL);   /* запомнить текущее время T2 */

    //delta_ms =  1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;

    printf("\nN=%d. Milliseconds passed: %f\n", N, delta_ms); /* T2 -T1 */
    printf("\nX: %f\n", X); /* T2 -T1 */
    return 0;
}
