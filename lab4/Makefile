CC=gcc
CFLAGS=-O3 -Wall -Werror

seq:
	$(CC) $(CFLAGS) -o lab4-seq -Wno-unknown-pragmas lab4.c -lm 
	
mp:
	$(CC) $(CFLAGS) -o lab4-mp -fopenmp lab4.c -lm 

	
test:
	$(CC) $(CFLAGS) -o test1 -fopenmp test.c 

