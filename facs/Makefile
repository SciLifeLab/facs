CFLAGS=-O3 -DFIFO -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE -Wall -fopenmp -g -DNODEBUG -lm -lz
.PHONY: clean 
.SUFFIXES:.c .o
PROG=facs
LOBJS= big_query.o bloom.o file_dir.o good_build.o lookup8.o suggestions.o tool.o simple_check_1_ge.o simple_remove.o simple_remove_l.o
AOBJS= big_query.o bloom.o file_dir.o good_build.o lookup8.o suggestions.o tool.o simple_check_1_ge.o simple_remove.o simple_remove_l.o

all:$(PROG)

mpi:
	@echo Make sure you have MPI support on your cluster hint: module load openmpi
	#mpicc -c *.c ${CFLAGS}
	#mpicc -c mpi_decon.c -O3 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE 
	#mpicc -o mpi_decon mpi_decon.o bloom.o suggestions.o lookup8.o  -lm ${CFLAGS}
	#mpicc -o mpi_bloom mpi_bloom.o bloom.o suggestions.o lookup8.o file_dir.o -lm ${CFLAGS}
	#mpicc -o mpi_bloom_l mpi_bloom_l.o bloom.o suggestions.o lookup8.o file_dir.o -lm ${CFLAGS}
	#mpirun -np 1 ./mpi_bloom_l -l tzcoolman  -q test.fna

.c.o:
		$(CC) -c $(DFLAGS) $(INCLUDES) $< -o $@ $(CFLAGS) 

clean:
	rm -f core.* vgcore.* *.o *.so *.a *.info ${PROG}

${PROG}:lib${PROG}.a $(AOBJS)

${PROG}:lib${PROG}.a $(AOBJS) main.o
		$(CC) $(DFLAGS) $(AOBJS) main.o -o $@ -L. -l${PROG} $(LIBS) $(CFLAGS) 

lib${PROG}.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

main.o: big_query.h bloom.h build.h check.h file_dir.h hashes.h tool.h remove.h remove_l.h