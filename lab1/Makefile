CC=gcc
CFLAGS=-O3 -Wall -Werror

seq:
	$(CC) $(CFLAGS) -o lab1-seq lab1.c -lm 

par:
	$(CC) $(CFLAGS) -floop-parallelize-all -ftree-parallelize-loops=2 lab1.c -o lab1-par-2 -lm
	$(CC) $(CFLAGS) -floop-parallelize-all -ftree-parallelize-loops=3 lab1.c -o lab1-par-3 -lm
	$(CC) $(CFLAGS) -floop-parallelize-all -ftree-parallelize-loops=4 lab1.c -o lab1-par-4 -lm
	$(CC) $(CFLAGS) -floop-parallelize-all -ftree-parallelize-loops=8 lab1.c -o lab1-par-8 -lm
icc_par:
	icc -fast -parallel -qopt-report-phase=par -par-threshold=2 -o lab1-icc-par-2 lab1.c
	icc -fast -parallel -qopt-report-phase=par -par-threshold=3 -o lab1-icc-par-3 lab1.c
	icc -fast -parallel -qopt-report-phase=par -par-threshold=4 -o lab1-icc-par-4 lab1.c
	icc -fast -parallel -qopt-report-phase=par -par-threshold=8 -o lab1-icc-par-8 lab1.c
	
icc_seq:
	icc -fast -o lab1-icc-seq lab1.c 
	
clang_seq:
	clang-10 -O3 -Wall -Werror -o lab1-seq-clang lab1.c -lm
	
clang_par:
	clang -O3 -Wall -Werror -mllvm -polly -mllvm -polly-parallel -lgomp -mllvm -polly-num-threads=2 -o lab1-par-clang-2 lab1.c -lm
	clang -O3 -Wall -Werror -mllvm -polly -mllvm -polly-parallel -lgomp -mllvm -polly-num-threads=3 -o lab1-par-clang-3 lab1.c -lm
	clang -O3 -Wall -Werror -mllvm -polly -mllvm -polly-parallel -lgomp -mllvm -polly-num-threads=4 -o lab1-par-clang-4 lab1.c -lm
	clang -O3 -Wall -Werror -mllvm -polly -mllvm -polly-parallel -lgomp -mllvm -polly-num-threads=8 -o lab1-par-clang-8 lab1.c -lm
