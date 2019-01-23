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
	rm -f sofos sofos_debug sofos_test sofos_test_coverage

SOFOSCC=sofos.cc main.cc
SOFOSHPP=sofos.hpp vcf.hpp

sofos: $(SOFOSCC) $(SOFOSHPP)
	$(CXX) $(CXXFLAGS) $(SOFOSFLAGS) -o $@ $(SOFOSCC)

sofos_debug: $(SOFOSCC) $(SOFOSHPP)
	$(CXX) -g -O0 $(SOFOSFLAGS) -o $@ $(SOFOSCC)

sofos_test: $(SOFOSCC) $(SOFOSHPP)
	$(CXX) -g -O0 $(SOFOSFLAGS) -DSOFOS_UNIT_TESTS -o $@ $(SOFOSCC)

sofos_test_coverage: $(SOFOSCC) $(SOFOSHPP)
	$(CXX) -g -O0 $(SOFOSFLAGS) $(GCOVFLAGS) -DSOFOS_UNIT_TESTS -o $@ $(SOFOSCC)

test: sofos_test
	./sofos_test

coverage: sofos_test_coverage
	rm -f sofos.gcda
	./sofos_test_coverage
	gcovr -r . -e catch.hpp --html-details -o coverage.html

.PHONY: coverage test
