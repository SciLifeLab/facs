CFLAGS=-O3 -DFIFO -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE -fopenmp -g -DDEBUG
.PHONY:  tests clean

all:
	@echo Make sure you have MPI support on your cluster hint: module load openmpi
	mpicc -c *.c ${CFLAGS}
	#mpicc -c mpi_decon.c -O3 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE 
	#mpicc -o mpi_bloom mpi_bloom.o bloom.o suggestions.o lookup8.o file_dir.o -lm ${CFLAGS}
	#mpicc -o mpi_bloom_l mpi_bloom_l.o bloom.o suggestions.o lookup8.o file_dir.o -lm ${CFLAGS}
	gcc -c simple_check_1_ge.c -O3
	${CC} -o raven main.o good_build.o simple_check_1_ge.o simple_remove.o simple_remove_l.o bloom.o suggestions.o lookup8.o file_dir.o tool.o -lm ${CFLAGS}
	${CC} -o test test.o good_build.o simple_check_1_ge.o simple_remove.o simple_remove_l.o bloom.o suggestions.o lookup8.o file_dir.o tool.o -lm ${CFLAGS}	
	#mpicc -o mpi_decon mpi_decon.o bloom.o suggestions.o lookup8.o  -lm ${CFLAGS}

classi:
	#mpicc -c *.c ${CFLAGS}
	gcc -c simple_remove_l.c ${CFLAGS}
	${CC} -o simple_remove_l simple_remove_l.o bloom.o suggestions.o file_dir.o lookup8.o -lm ${CFLAGS}

clean:
	rm -f core.* vgcore.* *.o *.so bloom_build simple_check simple_remove simple_remove_l

tests:
	mkdir -p tests/data
	test -s tests/data/ecoli_K12.fasta||wget http://togows.dbcls.jp/entry/ncbi-nucleotide/U00096.2.fasta -O tests/data/ecoli_K12.fasta
	#test -s tests/data/ecoli_dummy.fastq.gz||gzip tests/data/ecoli_dummy.fastq
	#./bloom_build -r tests/data/ecoli_K12.fasta -o tests/data/ecoli.bloom
	./tests/fastq_dummy.py 50 tests/data/ecoli_dummy.fastq
	test -s tests/data/ecoli_dummy.fastq.gz||gzip tests/data/ecoli_dummy.fastq
	test -e tests/data/ecoli_dummy.fastq.gz.fifo||mkfifo tests/data/ecoli_dummy.fastq.gz.fifo
	gunzip -c tests/data/ecoli_dummy.fastq.gz > tests/data/ecoli_dummy.fastq.gz.fifo &
	#./simple_check -m 1 -q tests/data/ecoli_dummy.fastq -f tests/data/ecoli.bloom
	./simple_check -m 1 -q tests/data/ecoli_dummy.fastq -f tests/data/ecoli_dummy.fastq.gz.fifo -r tests/data/ecoli.bloom

mpirun:
	mpirun -np 1 ./mpi_bloom_l -l tzcoolman  -q test.fna

classi_t:
	./simple_remove_l -q ~/test/test.fna -r tests/data/ecoli.bloom 

python:
	rm -rf build/ drass.so && python setup.py build_ext --inplace && python build.py
	gcc -c simple_check_1_ge.c -O3
