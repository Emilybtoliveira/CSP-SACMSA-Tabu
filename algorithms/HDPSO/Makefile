CPLEXFLAGS=-O3 -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD  -lm -pthread -std=c++0x -ldl

CFLAGS=-std=c++11 -static-libstdc++ -static-libgcc -Wall -g

all:HDPSO.o 
	g++  HDPSO.o -o HDPSO.run $(CPLEXFLAGS) $(CFLAGS)

%.o: %.cpp %.hpp
	g++ -c $< -o $@ $(CFLAGS)

HDPSO.o: HDPSO.cpp
	g++ -c -o HDPSO.o HDPSO.cpp $(CPLEXFLAGS) $(CFLAGS)

clean:
	rm -f *.o
