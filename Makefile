VERSION = 4.1
CC = gcc
ROOTSYS = /lib/root/root
CFLAGS = -Wall -g -I$(ROOTSYS)/include $(sh $ROOTSYS/bin/root-config --cflags)
LDFLAGS = $(sh $ROOTSYS/bin/root-config --cflags --glibs) -lgsl -lgslcblas -lm -lncurses -ggdb -lnest3 

OBJ = main.o


default: main

%.o: %.cpp
	$(CC) $(CFLAGS) -c $*.cpp

main: $(OBJ)
	$(CC)  $(CFLAGS) $(LDFLAGS) $(OBJ) 
