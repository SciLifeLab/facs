.PHONY: tests valgrind
PROG=facs

all:$(PROG)
	cd facs; make

tests: python
	nosetests -P -v -s tests

valgrind: python
	valgrind --tool=memcheck --suppressions=tests/utils/valgrind-python.supp nosetests -v -s tests/test_basic.py

python:
	rm -rf build/ ${PROG}.so && python setup.py build_ext && python setup.py install

clean:
	cd facs; make clean
