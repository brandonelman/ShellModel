CC = g++
CFLAGS = -g -Wall

all: ShellModel
ARMADILLO_BASE=/user/elman/libraries/armadillo-6.400.3
ARMADILLO_INCLUDE=$(ARMADILLO_BASE)/include

LIBS=-llapack -lblas -L$(ARMADILLO_BASE) 
INCLUDES=-Iinclude -I$(ARMADILLO_INCLUDE)          
#LIBS = -llapack -lblas -L/user/poxonpea/armadillo-6.500.5/ 

ShellModel:  lib/State.o lib/ShellModel.o

	$(CC) $(CFLAGS) $(INCLUDES) -o bin/shellmodel  lib/State.o lib/ShellModel.o $(LIBS)

lib/State.o: src/State.cpp include/State.hh
	$(CC) $(CFLAGS) $(INCLUDES) -c src/State.cpp -o lib/State.o $(LIBS)

lib/ShellModel.o: include/State.hh src/ShellModel.cpp include/ShellModel.hh
	$(CC) $(CFLAGS) $(INCLUDES) -c src/ShellModel.cpp -o lib/ShellModel.o $(LIBS)

clean: 
	$(RM) lib/*.o bin/*
