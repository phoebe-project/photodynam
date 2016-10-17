CXX=g++
#CXX=clang++
CXXFLAGS=-O3 -Wall -std=c++11 -I./include

all: photodynam

%: source/%.cpp
	$(CXX) $< -o $@ $(CXXFLAGS)
