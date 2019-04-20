CPP=g++
FLAGS=-Wall -O2


all:
	$(CPP) $(FLAGS) -o test benchmark.cpp

clean:
	rm test
