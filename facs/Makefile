# Common compiler flags for all OS
CFLAGS=-O3 -Wall -g -DNODEBUG

# Detect OS
UNAME := $(shell uname)
ifneq ($(UNAME), Darwin)
CFLAGS+=-fopenmp -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE
LIBS=-lm -lz
else
CFLAGS+=-Wno-unknown-pragmas -Wno-unused-value -arch x86_64
LIBS=-lz
endif

.PHONY: tests clean valgrind
.SUFFIXES:.c .o
PROG=facs
LOBJS= big_query.o bloom.o file_dir.o good_build.o lookup8.o suggestions.o tool.o simple_check_1_ge.o simple_remove.o prob.o info.o
AOBJS= big_query.o bloom.o file_dir.o good_build.o lookup8.o suggestions.o tool.o simple_check_1_ge.o simple_remove.o prob.o info.o

all:$(PROG)

tests: python
	nosetests -v -s ../tests

valgrind: python
	valgrind --tool=memcheck --suppressions=../tests/utils/valgrind-python.supp nosetests -v -s ../tests/test_basic.py

mpi:
	@echo Make sure you have MPI support on your cluster hint: module load openmpi
	#mpicc -c *.c ${CFLAGS} ${LIBS}
	#mpicc -c mpi_decon.c -O3 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE 
	#mpicc -o mpi_decon mpi_decon.o bloom.o suggestions.o lookup8.o  -lm ${CFLAGS}
	mpicc -c mpi_bloom.c ${CFLAGS}
	mpicc -o facs_mpi tool.o big_query.o simple_check_1_ge.o simple_remove.o mpi_bloom.o bloom.o suggestions.o lookup8.o file_dir.o prob.o ${LIBS} ${CFLAGS}

python:
	rm -rf build/ ${PROG}.so && python setup.py build_ext --inplace && python setup.py develop

clean:
	rm -f *.pyc core.* vgcore.* *.o *.so *.a *.info ${PROG}


.c.o:
		$(CC) -c $(DFLAGS) $(INCLUDES) $< -o $@ $(CFLAGS) 

${PROG}:lib${PROG}.a $(AOBJS)

${PROG}:lib${PROG}.a $(AOBJS) main.o
		$(CC) $(DFLAGS) $(AOBJS) main.o -o $@ -L. -l${PROG} $(LIBS) $(CFLAGS) 

lib${PROG}.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

main.o: big_query.h bloom.h build.h check.h file_dir.h hashes.h tool.h remove.h info.h

