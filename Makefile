CXX?=c++
CXXFLAGS?=-O2 -DNDEBUG
HTSLIB?=-lhts -lm
CLANGTIDY?=clang-tidy
CLANGFORMAT?=clang-format


SOFOSFLAGS=-std=c++11 $(HTSLIB)
GCOVFLAGS=-fprofile-arcs -ftest-coverage -fPIC
GPROFFLAGS=-pg

default: all

all: sofos

.PHONY: default all clean

clean:
	rm -f sofos sofos_debug sofos_test sofos_test_coverage
	rm -f *.gcov *.gcda *.gcno
	rm -f coverage*.html

SOFOSCC=sofos.cc main.cc
SOFOSHPP=sofos.hpp vcf.hpp
FORMATFILES=$(SOFOSCC) $(SOFOSHPP) test_utils.hpp

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
	rm -f *.gcda
	./sofos_test_coverage
	gcovr -r . -e catch.hpp -e test_utils.hpp --html-details -o coverage.html

test_codecov: sofos_test_coverage
	rm -f *.gcda
	./sofos_test_coverage
	for filename in $(SOFOSCC); do gcov -n -o . $$filename > /dev/null; done
	bash -c 'bash <(curl -s https://codecov.io/bash)'

tidy:
	$(CLANGTIDY) $(SOFOSCC) -- $(CXXFLAGS) $(SOFOSFLAGS) -Wall

format:
	$(CLANGFORMAT) -i $(SOFOSCC) $(SOFOSHPP) $(SOFOSTESTHPP)

test_format: sofos_test_coverage
	test 0 -eq `$(CLANGFORMAT) -output-replacements-xml $(FORMATFILES) 2>/dev/null | grep offset | wc -l`

test_tidy: sofos_test_coverage
	test 0 -eq `$(CLANGTIDY) $(SOFOSCC) -- $(CXXFLAGS) $(SOFOSFLAGS) -Wall 2>/dev/null | wc -l`

.PHONY: coverage tidy format test test_codecov test_format
