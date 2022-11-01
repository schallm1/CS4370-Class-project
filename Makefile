# Makefile for clothing distribution base code
# Michael Schall


CC      = g++
CFLAGS  = -c
TARGET  = P1
OBJS    = P1.o

$(TARGET):	$(OBJS)
		$(CC) -o $(TARGET) $(OBJS)

P1.o:		P1.cpp
		$(CC) $(CFLAGS)  P1.cpp -std=c++11

clean:
		/bin/rm -f *.o $(TARGET)
