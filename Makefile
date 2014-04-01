.PHONY: tests valgrind
PROG=facs

all:$(PROG)
	cd facs; make
mpi: 
	cd facs; make; make mpi
tests: python
	cd tests && nosetests -v -s -P --with-timer

benchmarks: python
	# XXX: Add a helper function that determines the number of cores present
	cd tests &&  \
	nosetests -v -s -P --with-timer test_simngs.py && \
	export OMP_NUM_THREADS=1 && nosetests -v -s -P --with-timer test_basic.py && \
	export OMP_NUM_THREADS=8 && nosetests -v -s -P --with-timer test_basic.py && \
	export OMP_NUM_THREADS=16 && nosetests -v -s -P --with-timer test_basic.py && \
	export OMP_NUM_THREADS=1 && nosetests -v -s -P --with-timer test_fastqscreen.py && \
	export OMP_NUM_THREADS=8 && nosetests -v -s -P --with-timer test_fastqscreen.py && \
	export OMP_NUM_THREADS=16 && nosetests -v -s -P --with-timer test_fastqscreen.py

valgrind: python
	valgrind --tool=memcheck --suppressions=facs/utils/valgrind-python.supp nosetests -P -v -s

python: clean
	ARCHFLAGS=-Wno-error=unused-command-line-argument-hard-error-in-future python setup.py install

clean:
	rm -rf build dist $(PROG).egg-info
	cd facs; make clean
