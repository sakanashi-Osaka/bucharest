target=decay_Dron.cpp
object=decay_Dron.o
output=decay_Dron
#target=decay13C.cpp
#object=decay13C.o
#output=decay13C

CC = g++ 

ROOT = $(shell root-config --cflags --libs)

all: main

main: decay_Dron.o lorlib.o lorlib.h kinelib.o kinema.h
	$(CC) $(ROOT) $(object) -o $(output) 

decay_Dron.o: $(target)
	$(CC) $(ROOT) -c $(target)

lorlib.o: lorlib.cpp lorlib.h
	$(CC) -c lorlib.cpp

kinelib.o: kinelib.cpp kinelib.h kinema.h
	$(CC) -c kinelib.cpp

clean:
	rm -f decay_Dron *.o
