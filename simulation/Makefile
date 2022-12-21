target=decay.cpp
output=decay

CC = g++ 

ROOT = $(shell root-config --cflags --libs)

all: main

main: decay.o lorlib.o lorlib.h kinelib.o kinema.h
	$(CC) $(ROOT) decay.o -o $(output) 

decay.o: $(target)
	$(CC) $(ROOT) -c $(target)

lorlib.o: lorlib.cpp lorlib.h
	$(CC) -c lorlib.cpp

kinelib.o: kinelib.cpp kinelib.h kinema.h
	$(CC) -c kinelib.cpp

clean:
	rm -f decay *.o
