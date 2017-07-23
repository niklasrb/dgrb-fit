VERSION = 4.1
CC = g++
#ROOTSYS = /lib/root/root

CFLAGS = -Wall -g  $(sh $ROOTSYS/bin/root-config --cflags) #-I$(ROOTSYS)/include
LDFLAGS = -ggdb $(sh $ROOTSYS/bin/root-config --cflags --glibs) 
#-lgsl -lgslcblas -lncurses -L$(ROOTSYS)/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic

#$(sh $ROOTSYS/bin/root-config --cflags --glibs) -lgsl -lgslcblas -lm -lncurses -ggdb #-lnest3 

SOURCE = main.cpp
OBJ = main.o


default: main

%.o: %.cpp
	$(CC) $(CFLAGS) -c $*.cpp

main: $(OBJ)
	$(CC)  $(CFLAGS) $(LDFLAGS) $(SOURCE) 
