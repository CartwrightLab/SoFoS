CXX?=c++
DEBUG?=0
HTSLIB?=-lhts -lm

ifeq ($(DEBUG), 1)
	CXXFLAGS?=-g -O0
else
	CXXFLAGS?=-O2 -DNDEBUG
endif
CXXFLAGS+=-std=c++11

default: all

all: sofos

.PHONY: default all clean

clean:
	rm -f sofos

sofos: sofos.cc
	$(CXX) $(CXXFLAGS) $(HTSLIB) -o $@ $<

