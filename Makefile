CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE -fopenmp -g
.PHONY: tests clean

all:
	gcc -c *.c ${CFLAGS}
	gcc -o bloom_build good_build.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}
	gcc -o simple_check simple_check_1_ge.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}
	gcc -o simple_remove simple_remove.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}

clean:
	rm *.o

tests:
	./tests/fastq_dummy.py 50 tests/data/ecoli_dummy.fastq
	./bloom_build -k 21 -l tests/reference_genomes.list -p tests/data/ecoli.bloom
	./simple_check -m 1 -q tests/data/ecoli_dummy.fastq -l tests/bloom_filters.list
