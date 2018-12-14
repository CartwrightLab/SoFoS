CXX:=c++

default: all

all: sofos

.PHONY: default all

sofos: sofos.cc
	$(CXX) -O2 -std=c++11 -lhts -lm -o $@ $<