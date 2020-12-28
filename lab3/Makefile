CC=gcc
CFLAGS=-O3 -Wall -Werror

seq:
	$(CC) $(CFLAGS) -o lab3-seq -Wno-unknown-pragmas lab3.c -lm 
	
mp: lab3.c
	$(CC) $(CFLAGS) -o lab3-mp -fopenmp lab3.c -lm 
	
sch:
	$(CC) $(CFLAGS) -o lab3-mp-schedule -fopenmp lab3_schedule.c -lm 
	
	


