CXX?=c++
CXXFLAGS?=-O2 -DNDEBUG
CLANGTIDY?=clang-tidy
CLANGFORMAT?=clang-format

LOCAL_HTSLIB_VERSION=1.9

ifdef USE_LOCAL_HTSLIB
HTSLIBFLAGS?=-I./htslib-$(LOCAL_HTSLIB_VERSION)
HTSLIBS?=./htslib-$(LOCAL_HTSLIB_VERSION)/libhts.a -lz -lm -llzma -lbz2 -pthread
else
HTSLIBS?=-lhts -lm
endif

SOFOSFLAGS=-std=c++11 $(HTSLIBFLAGS)
SOFOSLIBS=$(HTSLIBS)
GCOVFLAGS=-fprofile-arcs -ftest-coverage -fPIC
GPROFFLAGS=-pg

default: all

all: sofos

.PHONY: default all clean

clean:
	rm -f sofos sofos_debug unittest unittest_coverage sofos_coverage
	rm -f *.o
	rm -f *.gcov *.gcda *.gcno
	rm -f coverage*.html cobertura.xml
	rm -f clang-tidy.txt clang-format.patch

SOFOSCC=sofos.cc main.cc
UNITCC=sofos.cc unittest.cc
SOFOSHPP=sofos.hpp vcf.hpp
FORMATFILES=$(SOFOSCC) $(SOFOSHPP) unittest.cc
TIDYFILES=$(SOFOSCC) unittest.cc

sofos: $(SOFOSCC) $(SOFOSHPP)
	$(CXX) $(CXXFLAGS) $(SOFOSFLAGS) -o $@ $(SOFOSCC) $(SOFOSLIBS)

sofos_debug: $(SOFOSCC) $(SOFOSHPP)
	$(CXX) -g -O0 $(SOFOSFLAGS) -o $@ $(SOFOSCC) $(SOFOSLIBS)

unittest: $(UNITCC) $(SOFOSHPP)
	$(CXX) -g -O0 $(SOFOSFLAGS) -DSOFOS_UNIT_TESTS -o $@ $(UNITCC) $(SOFOSLIBS)

%_cov.o: %.cc
	$(CXX) -c -g -O0 $(SOFOSFLAGS) $(GCOVFLAGS) -o $@ $<

sofos_cov.o : $(SOFOSHPP)

unit_cov.o : $(SOFOSHPP)

sofos_coverage: $(patsubst %.cc,%_cov.o,$(SOFOSCC))
	$(CXX) -g -O0 $(SOFOSFLAGS) $(GCOVFLAGS) -o $@ $^ $(SOFOSLIBS)

unittest_coverage: $(patsubst %.cc,%_cov.o,$(UNITCC))
	$(CXX) -g -O0 $(SOFOSFLAGS) $(GCOVFLAGS) -DSOFOS_UNIT_TESTS -o $@ $^ $(SOFOSLIBS)

test: unittest sofos_debug
	./unittest
	cd test && bash run_tests.bash ../sofos_debug

coverage: sofos_coverage unittest_coverage
	rm -f *.gcda
	./unittest_coverage
	cd test && bash run_tests.bash ../sofos_coverage	

coverage_html: coverage
	gcovr -r . -s -e contrib/catch.hpp -e unittest.cc --html-details -o coverage.html

coverage_xml: coverage
	gcovr -r . -s -e contrib/catch.hpp -e unittest.cc -x -o cobertura.xml

codecov.bash:
	curl -s https://codecov.io/bash > $@

# edit the xml inplace to upgrade partial hits to full hits
test_codecov: coverage_xml codecov.bash
	sed -i -e 's|condition-coverage="[^0"][^"]*\([0-9]\+\))"|condition-coverage="100% (\1/\1)"|g' cobertura.xml
	bash codecov.bash -X gcov search

tidy:
	$(CLANGTIDY) $(TIDYFILES) -- $(CXXFLAGS) $(SOFOSFLAGS) -Wall

format:
	$(CLANGFORMAT) -i $(SOFOSCC) $(SOFOSHPP) $(FORMATFILES)

test_format:
	@echo 'Testing for code format issues...'
	@bash -c 'diff -u <(cat $(FORMATFILES)) <($(CLANGFORMAT) $(FORMATFILES) 2>/dev/null) | tee clang-format.patch'
	@test 0 -eq `cat clang-format.patch | wc -l`

test_tidy:
	@echo 'Analyzing code for potential problems...'
	@$(CLANGTIDY) $(SOFOSCC) -- $(CXXFLAGS) $(SOFOSFLAGS) -Wall 2>/dev/null | tee clang-tidy.txt
	@test 0 -eq `cat clang-tidy.txt | wc -l`

.PHONY: coverage tidy format test test_codecov test_format test_tidy coverage_html coverage_xml


htslib-$(LOCAL_HTSLIB_VERSION).tar.bz2:
	curl -SL https://github.com/samtools/htslib/releases/download/$(LOCAL_HTSLIB_VERSION)/htslib-$(LOCAL_HTSLIB_VERSION).tar.bz2 > $@

htslib-$(LOCAL_HTSLIB_VERSION): htslib-$(LOCAL_HTSLIB_VERSION).tar.bz2
	tar xf $<

htslib: htslib-$(LOCAL_HTSLIB_VERSION)
	cd htslib-$(LOCAL_HTSLIB_VERSION) && $(MAKE)

.PHONY: htslib
