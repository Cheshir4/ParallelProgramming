all: build

build: *.c
	gcc -m64 -Wall -Werror -o lab2-fw $< -I ~/Загрузки/FW_1.3.1_Lin64 -L ~/Загрузки/FW_1.3.1_Lin64/lib/ -lm -lfwImage -lfwBase -lfwSignal

run: build
	LD_LIBRARY_PATH=~/Загрузки/FW_1.3.1_Lin64/lib ./lab2-fw 1064 1