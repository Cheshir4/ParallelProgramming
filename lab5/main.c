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
        struct lab_data* data = (struct lab_data*)arg;
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
    return NULL;
}

int main(int argc, char *argv[]) {

    if ( pthread_mutex_init( &mutex, NULL) != 0 )
        printf( "mutex init failed\n" );

    struct lab_data labData;

    float X = 0;

    struct timeval T1, T2;
    long delta_ms, delta_ms_Generate=0, delta_ms_Map=0, delta_ms_Merge=0, delta_ms_Sort=0, delta_ms_Reduce=0;
    long time_stamp_Generate=0, time_stamp_Map=0, time_stamp_Merge=0, time_stamp_Sort=0, time_stamp_Reduce=0;

    gettimeofday(&T1, NULL); /* запомнить текущее время T1 */
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

        gettimeofday(&T2, NULL);   /* запомнить текущее время T2 */
        time_stamp_Generate =  1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
        delta_ms_Generate += (time_stamp_Generate - time_stamp_Reduce);

        PARALLEL_START(map_m1, labData, 2, labData.N)
        PARALLEL_START(map_m2, labData, 2, labData.N / 2)

        gettimeofday(&T2, NULL);   /* запомнить текущее время T2 */
        time_stamp_Map =  1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
        delta_ms_Map += (time_stamp_Map - time_stamp_Generate);

        PARALLEL_START(merge, labData, 2, labData.N / 2)

        gettimeofday(&T2, NULL);   /* запомнить текущее время T2 */
        time_stamp_Merge =  1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
        delta_ms_Merge += (time_stamp_Merge - time_stamp_Map);


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

        gettimeofday(&T2, NULL);   /* запомнить текущее время T2 */
        time_stamp_Sort =  1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
        delta_ms_Sort += (time_stamp_Sort - time_stamp_Merge);


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

        gettimeofday(&T2, NULL);   /* запомнить текущее время T2 */
        time_stamp_Reduce =  1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
        delta_ms_Reduce += (time_stamp_Reduce - time_stamp_Sort);

        free(labData.m1);
        free(labData.m2);
        free(labData.m2_copy);
    }


    gettimeofday(&T2, NULL);   /* запомнить текущее время T2 */
    delta_ms =  1000*(T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec)/1000;
    printf("\nN=%d. Milliseconds passed after Generate: %ld\n", labData.N, delta_ms_Generate); /* T2 -T1 */
    printf("\nN=%d. Milliseconds passed after Map: %ld\n", labData.N, delta_ms_Map); /* T2 -T1 */
    printf("\nN=%d. Milliseconds passed after Merge: %ld\n", labData.N, delta_ms_Merge); /* T2 -T1 */
    printf("\nN=%d. Milliseconds passed after Sort: %ld\n", labData.N, delta_ms_Sort); /* T2 -T1 */
    printf("\nN=%d. Milliseconds passed after Reduce: %ld\n", labData.N, delta_ms_Reduce); /* T2 -T1 */
    printf("\nN=%d. Milliseconds passed: %ld\n", labData.N, delta_ms); /* T2 -T1 */
    pthread_join(percent_thread, NULL);

    printf("X = %.6f\n", X);
    pthread_mutex_destroy(&mutex);


    return 0;
}
