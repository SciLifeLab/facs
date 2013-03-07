.PHONY: tests valgrind
PROG=facs

all:$(PROG)
	cd facs; make

tests: python
	nosetests -P -v -s tests

valgrind: python
	valgrind --tool=memcheck --suppressions=facs/utils/valgrind-python.supp nosetests -P -v -s

python:
	rm -rf build/ ${PROG}.so && python setup.py build_ext && python setup.py install

clean:
	cd facs; make clean
