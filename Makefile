.PHONY: tests valgrind
PROG=facs

all:$(PROG)
	cd facs; make

tests: python
	#XXX add standard suite, do not run long running jobs on Travis-CI
	cd tests && nosetests -v -s -P --with-timer

valgrind: python
	valgrind --tool=memcheck --suppressions=facs/utils/valgrind-python.supp nosetests -P -v -s

python: clean
	python setup.py install

clean:
	rm -rf build dist $(PROG).egg-info
	cd facs; make clean
