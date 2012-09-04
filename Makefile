CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE -fopenmp# -g -DDEBUG
.PHONY: tests clean

all:
	${CC} -c *.c ${CFLAGS}
	${CC} -o bloom_build good_build.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}
	${CC} -o simple_check simple_check_1_ge.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}
	${CC} -o simple_remove simple_remove.o bloom.o suggestions.o lookup8.o -lm ${CFLAGS}

clean:
	rm -f core.* *.o bloom_build simple_check simple_remove

tests:
	mkdir -p tests/data
	test -s tests/data/ecoli_K12.fasta||wget http://togows.dbcls.jp/entry/ncbi-nucleotide/U00096.2.fasta -O tests/data/ecoli_K12.fasta
	./tests/fastq_dummy.py 50 tests/data/ecoli_dummy.fastq
	test -s tests/data/ecoli_dummy.fastq.gz||gzip tests/data/ecoli_dummy.fastq
	./bloom_build -r tests/data/ecoli_K12.fasta -o tests/data/ecoli.bloom
	test -e tests/data/ecoli_dummy.fastq.gz.fifo||mkfifo tests/data/ecoli_dummy.fastq.gz.fifo
	gunzip -c tests/data/ecoli_dummy.fastq.gz > tests/data/ecoli_dummy.fastq.gz.fifo &
	./simple_check -m 1 -q tests/data/ecoli_dummy.fastq -r tests/data/ecoli.bloom
	@echo Checking contamination against gz fifo file...
	./simple_check -m 1 -q tests/data/ecoli_dummy.fastq.gz.fifo -r tests/data/ecoli.bloom
