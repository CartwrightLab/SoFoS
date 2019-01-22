CXX?=c++
CXXFLAGS?=-O2 -DNDEBUG
HTSLIB?=-lhts -lm

SOFOSFLAGS=-std=c++11 $(HTSLIB)
GCOVFLAGS=-fprofile-arcs -ftest-coverage -fPIC
GPROFFLAGS=-pg

default: all

all: sofos

.PHONY: default all clean

clean:
	rm -f sofos

sofos: sofos.cc
	$(CXX) $(CXXFLAGS) $(SOFOSFLAGS) -o $@ $^

sofos_debug: sofos.cc
	$(CXX) -g -O0 $(SOFOSFLAGS) -o $@ $^

sofos_test: sofos.cc
	$(CXX) -g -O0 $(SOFOSFLAGS) -DSOFOS_UNIT_TESTS -o $@ $^

sofos_test_coverage: sofos.cc
	$(CXX) -g -O0 $(SOFOSFLAGS) $(GCOVFLAGS) -DSOFOS_UNIT_TESTS -o $@ $^

test: sofos_test
	./sofos_test

coverage: sofos_test_coverage
	rm -f sofos.gcda
	./sofos_test_coverage
	gcovr -r . -e catch.hpp --html-details -o coverage.html

.PHONY: coverage test
