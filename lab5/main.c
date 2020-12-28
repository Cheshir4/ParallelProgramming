#include "main.h"

#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>

#define A 480

const float e = 0.00001;

int current = 0;
int step = 2;
int length = 0;

pthread_mutex_t mutex;

int add(int local) {
    int result;
    pthread_mutex_lock(&mutex);

    result = current;
    current += local;

    pthread_mutex_unlock(&mutex);
    return result;
}

PARALLEL(generate_m1, {
    data->m1[i] = rand_r(data->my_seed) % A + 1;
})

PARALLEL(generate_m2, {
    data->m2[i] = rand_r(data->my_seed) % (9 * A + 1) + A;
    //printf("m2[%d]: %f \n", j, m2[j]);
    /* Подготавливаем копию для удобства сложения */
    data->m2_copy[i + 1] = data->m2[i];
})

PARALLEL(map_m1, {
    data->m1[i] = 1 / tanh(sqrt(data->m1[i]));
})

PARALLEL(map_m2, {
    if (i != 0) {
        data->m2[i] = data->m2[i] + data->m2_copy[i];
    }

    data->m2[i] = fabs(tan(data->m2[i]));
})

PARALLEL(merge, {
    data->m2[i] = (data->m1[i] < data->m2[i]) ? data->m1[i] : data->m2[i];
})

PARALLEL(sort, {
    int j = start;
    //printf("%i\n", i);
    while (j < end) {
        if (data->m2[j + 1] < data->m2[j]) {
            int tmp = data->m2[j];
            data->m2[j] = data->m2[j + 1];
            data->m2[j + 1] = tmp;
            j = i;
        } else j++;
    }
})

void* parallel_reduce(void* arg) {
        struct lab_data* data = (struct labData*)arg;
        int start = add(step);
        int end = start + step;
        if(end > length) end = length;
        float* local_x = (float*)malloc(sizeof(float));
        *local_x = 0;
        while (start < length) {
            for(int i = start; i < end; i ++) {
                int temp = 0;

                temp = (int) (data->m2[i] / data->min);

                if (temp % 2 == 0) {
                    *local_x += sin(data->m2[i]);
                }
            }
            start = add(step);
            end = start + step;
            if(end > length) end = length;
        }
        return (void*)local_x;
}

int it = 0;
void* percent(void* arg) {
    printf("Progress: %3d%%\n", it * 2);
    while(it < 50) {
        sleep(1);
        printf("Progress: %3d%%\n", it * 2);
    }
}

int main(int argc, char *argv[]) {

    if ( pthread_mutex_init( &mutex, NULL) != 0 )
        printf( "mutex init failed\n" );

    struct lab_data labData;

    float X = 0;

    pthread_t percent_thread;
    pthread_create(&percent_thread, NULL, percent, NULL);

    for(it = 0; it < 50; it++) {

        labData.my_seed[0] = it;

        labData.N = atoi(argv[1]); /* N равен первому параметру командной строки */
        labData.M = atoi(argv[2]) + 1;

        labData.m1 = (float *) malloc(sizeof(float) * labData.N);
        labData.m2 = (float *) malloc(sizeof(float) * labData.N / 2);
        labData.m2_copy = (float *) malloc(sizeof(float) * labData.N / 2 + 1);
        labData.m2_copy[0] = 0;

        PARALLEL_START(generate_m1, labData, 2, labData.N)
        PARALLEL_START(generate_m2, labData, 2, labData.N / 2)
        PARALLEL_START(map_m1, labData, 2, labData.N)
        PARALLEL_START(map_m2, labData, 2, labData.N / 2)
        PARALLEL_START(merge, labData, 2, labData.N / 2)
        PARALLEL_START(sort, labData, 2, labData.N / 2 - 2)

        int j = 0;
        while (j < (labData.N / 2) - 1) {
            if (labData.m2[j + 1] < labData.m2[j]) {
                int tmp = labData.m2[j];
                labData.m2[j] = labData.m2[j + 1];
                labData.m2[j + 1] = tmp;
                j = 0;
            } else j++;
        }

        labData.min = 0;

        while (labData.min == 0) {
            labData.min = (fabs(labData.m2[j]) < e) ? labData.m2[j] : 0;
            j++;
        }

        {
            int threads_num = labData.M;
            pthread_t threads[threads_num];
            current = 0;
            step = 2;
            length = labData.N / 2;
            for (int i = 1; i < threads_num; i++) {
                pthread_create(&threads[i], NULL, parallel_reduce, &labData);
            }
            for (int i = 1; i < threads_num; i++) {
                void *local_x;
                pthread_join(threads[i], &local_x);
                X += *((float *) local_x);
                free(local_x);
            }
        }

        free(labData.m1);
        free(labData.m2);
        free(labData.m2_copy);
    }

    pthread_join(percent_thread, NULL);

    printf("X = %.6f", X);
    pthread_mutex_destroy(&mutex);

    return 0;
}
