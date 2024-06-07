CXX=g++
#CXX=clang++
CXXFLAGS=-O3 -Wall -std=c++11 -I./include

all: photodynam

%: photodynam/%.cpp
	$(CXX) $< -o $@ $(CXXFLAGS)
