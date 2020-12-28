//
// Created by cheshir on 28.12.2020.
//

#ifndef LAB5_MAIN_H
#define LAB5_MAIN_H

struct lab_data {
    int N;
    int M;
    unsigned int my_seed[1];

    float *m1;
    float *m2;
    float *m2_copy;

    float min;
};

#define PARALLEL(name, work) \
void* parallel_##name(void* arg) { \
        struct lab_data* data = (struct labData*)arg; \
        int start = add(step); \
        int end = start + step; \
        if(end > length) end = length;                \
        while (start < length) {   \
            for(int i = start; i < end; i ++) \
                work \
            start = add(step);\
            end = start + step;\
            if(end > length) end = length;\
        }\
        return NULL;                     \
}

#define PARALLEL_START(name, work, step_val, size_val) \
{ \
    int threads_num = work.M; \
    pthread_t threads[threads_num]; \
    current = 0; \
    step = step_val; \
    length = size_val; \
    for (int i = 1; i < threads_num; i++) { \
        pthread_create(&threads[i], NULL, parallel_##name, &work); \
    } \
    for (int i = 1; i < threads_num; i++) { \
        pthread_join(threads[i], NULL); \
    } \
}

#endif //LAB5_MAIN_H
