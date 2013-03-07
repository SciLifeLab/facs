.PHONY: tests valgrind
PROG=facs

all:$(PROG)
	cd facs; make

tests: python
	cd tests && nosetests -v -s -P

valgrind: python
	valgrind --tool=memcheck --suppressions=facs/utils/valgrind-python.supp nosetests -P -v -s

python: clean
	python setup.py build_ext && python setup.py install

clean:
	rm -rf build dist $(PROG).egg-info
	cd facs; make clean
