CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE -fopenmp

all:
	gcc -c *.c ${CFLAGS}
	gcc -o bloom_build good_build.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}
	gcc -o simple_check simple_check_1_ge.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}
	gcc -o simple_remove simple_remove.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}

clean:
	rm *.o
