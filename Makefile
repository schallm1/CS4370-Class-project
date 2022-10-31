
all: output

output: P1.o
	g++ -o output P1.o

P1.o: P1.cpp
	g++ -c P1.cpp

clean:
	rm -f P1.o
	rm -f output
