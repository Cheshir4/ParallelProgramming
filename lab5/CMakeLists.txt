cmake_minimum_required(VERSION 3.17)
project(lab5 C)

set(CMAKE_C_STANDARD 11)

add_executable(lab5 main.c main.h)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(lab5 pthread)
target_link_libraries(lab5 m)