.PHONY: tests valgrind
PROG=facs

all:$(PROG)
	cd facs; make
mpi: 
	cd facs; make; make mpi
tests: python
	cd tests && nosetests -v -s -P --with-timer

benchmarks: python
	cd tests &&  \
	nosetests -v -s --with-timer test_basic.py && \
	nosetests -v -s --with-timer test_simngs.py && \
	nosetests -v -s --with-timer test_fastqscreen.py && \
	nosetests -v -s --with-timer test_deconseq.py

valgrind: python
	valgrind --tool=memcheck --suppressions=facs/utils/valgrind-python.supp nosetests -P -v -s

python: clean
	python setup.py install

clean:
	rm -rf build dist $(PROG).egg-info
	cd facs; make clean
