CPLEXFLAGS=-O3 -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD  -lm -pthread -std=c++0x -ldl

CFLAGS=-std=c++11 -static-libstdc++ -static-libgcc -Wall -g

all:SAGW.o 
	g++  SAGW.o -o SAGW.run $(CPLEXFLAGS) $(CFLAGS)

%.o: %.cpp %.hpp
	g++ -c $< -o $@ $(CFLAGS)

SAGW.o: SAGW.cpp
	g++ -c -o SAGW.o SAGW.cpp $(CPLEXFLAGS) $(CFLAGS)

clean:
	rm -f *.o
