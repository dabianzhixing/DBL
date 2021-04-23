CXX = g++
CXXFLAGS = -g -Wall -Wextra -fopenmp -O3 -std=c++11 -I /home/lqy/cuda-workspace/boost_1_73_0 -L /home/lqy/cuda-workspace/boost_1_73_0/stage/lib

all: construct query insert 


construct: construct.cc DBL.cc DBL.h
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_serialization

query: query.cc DBL.cc DBL.h
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_serialization
	
insert: insert.cc DBL.cc DBL.h
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_serialization

