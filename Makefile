CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE -fopenmp
.PHONY: tests clean

# If on linux:
#CC=gcc

# If on MacOS, make sure that you have a recent enough compiler that
# has support for OpenMP:
CC=gcc-mp-4.7


all:
	${CC} -c *.c ${CFLAGS}
	${CC} -o bloom_build good_build.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}
	${CC} -o simple_check simple_check_1_ge.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}
	${CC} -o simple_remove simple_remove.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}

clean:
	rm *.o

tests:
	mkdir -p tests/data
	test -s tests/data/ecoli_K12.fasta||wget http://togows.dbcls.jp/entry/ncbi-nucleotide/eschColi_K12,U00096.2.fasta -O tests/data/ecoli_K12.fasta
	./tests/fastq_dummy.py 50 tests/data/ecoli_dummy.fastq
	./bloom_build -k 21 -l tests/reference_genomes.list -p tests/data/ecoli.bloom
	./simple_check -m 1 -q tests/data/ecoli_dummy.fastq -l tests/bloom_filters.list
