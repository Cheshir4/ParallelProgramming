#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define A 480

int main(int argc, char* argv[])
{
    int i, N, j, tmp;
    float X = 0;
    struct timeval T1, T2;
    long delta_ms;
    float e = 0.00001;

    N = atoi(argv[1]); /* N равен первому параметру командной строки */
    gettimeofday(&T1, NULL); /* запомнить текущее время T1 */
    for (i=0; i<50; i++) /* 50 экспериментов */
        {

            srand(i);  /* инициализировать начальное значение ГСЧ   */
            //unsigned int *my_seed = malloc(sizeof(unsigned int));
            unsigned int my_seed[1];
            my_seed[0] = i;
            /* Заполнить массив исходных данных размером N */
            /* Решить поставленную задачу, заполнить массив с результатами*/
            float m1[N];
            for (j=0; j<N; j++)
                {
                    m1[j] = rand_r(my_seed) % A + 1; 
                    m1[j] = 1 / tanh(sqrt(m1[j]));
                }
                
            float m2[N/2];
            float m2_copy[N/2];
            for (j=0; j<N/2; j++)
                {
                    m2[j] = rand_r(my_seed) % (9*A + 1) + A;
                    m2_copy[j] = m2[j];
                    if (j != 0)
                        {
                           m2[j] = m2[j] + m2_copy[j-1];
                        } 

                    m2[j] = fabs(tan(m2[j]));
                }


            /* Этап Merge */
            for (j=0; j<N/2; j++)
                {
                    m2[j] = (m1[j] < m2[j]) ? m1[j] : m2[j];
                }

            /* Отсортировать массив с результатами указанным методом */

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
            /* Этап Reduce */

            float min = 0;
            j = 0;
            tmp = 0;

            while (min == 0)
                {
                    min = (fabs(m2[j]) < e) ? m2[j] : 0;
                    j++;
                }

            for (j=0; j<N/2; j++)
                {
                    tmp = (int)(m2[j] / min);

                    if (tmp % 2 == 0)
                        {
                            X = X + sin(m2[j]);
                        }
                }

        }
        gettimeofday(&T2, NULL);   /* запомнить текущее время T2 */
        delta_ms =  1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
        printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms); /* T2 -T1 */
        printf("\nX: %f\n", X); /* T2 -T1 */
        return 0;
}
