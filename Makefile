VERSION := 4.1
CC := g++

DEBUG := TRUE

#ROOTSYS = /usr/lib/root
LIBNEST := /usr/local/lib/libnest/

CFLAGS := -Wall -g -I$(ROOTSYS)/include $(shell $(ROOTSYS)/bin/root-config --cflags) 
LDFLAGS := $(shell $(ROOTSYS)/bin/root-config --cflags --glibs) -lgsl -lgslcblas -L$(LIBNEST) -lnest3 -llapack -lgfortran

ifeq ($(DEBUG), TRUE)
	LDFLAGS := -v -ggdb $(LDFLAGS) 
endif

DEPS :=  Constants.h InterpolationWrapper.h  CosmologyModel.h Benchmarking.h HaloModel.h Source.h DarkMatter.h AstrophysicalSource.h BLLAC.h FSRQ.h MAGN.h SFG.h EBLAbsorbtionCoefficient.h GalaxyCatalog.h AngularPowerSpectrum.h LinearMatterPowerSpectrum.h LoadFromFiles.h
SOURCE := main.cpp


default: main


main: $(DEPS) $(SOURCE)
	$(CC) -o main -O3 $(SOURCE) $(CFLAGS) $(LDFLAGS) 

clear:
	rm main 
