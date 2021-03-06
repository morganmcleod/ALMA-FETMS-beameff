# Project: beameff_64
# Makefile created by Dev-C++ 4.9.9.2
# then modified by T. Hunter for linux
# 1.3.6  MTM using iniparser 4.0 (from ../iniparser) instead of prior 1.8. 

CPP  = g++
CC   = gcc -DLINUX -g
OBJ  = dictionary.o iniparser.o brent.o f1dim.o mnbrak.o constants.o linmin.o frprmn.o nrutil.o fitamplitude.o nsi.o pointingangles.o plotcircles.o plotting_copol.o plotting_crosspol.o main.o getarrays.o fitphase.o efficiency.o outputfilefunctions.o SCANDATA.o utilities.o z.o $(RES)
LIBS = -lm
INC = -I../iniparser/src
CFLAGS = $(INC)
BIN  = beameff_64
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before beameff_64 all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o "beameff_64" $(LIBS)

dictionary.o: ../iniparser/src/dictionary.c
	$(CC) -c ../iniparser/src/dictionary.c -o dictionary.o $(CFLAGS)

iniparser.o: ../iniparser/src/iniparser.c
	$(CC) -c ../iniparser/src/iniparser.c -o iniparser.o $(CFLAGS)

brent.o: brent.c
	$(CC) -c brent.c -o brent.o $(CFLAGS)

f1dim.o: f1dim.c
	$(CC) -c f1dim.c -o f1dim.o $(CFLAGS)

mnbrak.o: mnbrak.c
	$(CC) -c mnbrak.c -o mnbrak.o $(CFLAGS)

constants.o: constants.c
	$(CC) -c constants.c -o constants.o $(CFLAGS)

linmin.o: linmin.c
	$(CC) -c linmin.c -o linmin.o $(CFLAGS)

frprmn.o: frprmn.c
	$(CC) -c frprmn.c -o frprmn.o $(CFLAGS)

nrutil.o: nrutil.c
	$(CC) -c nrutil.c -o nrutil.o $(CFLAGS)

fitamplitude.o: fitamplitude.c
	$(CC) -c fitamplitude.c -o fitamplitude.o $(CFLAGS)

nsi.o: nsi.c
	$(CC) -c nsi.c -o nsi.o $(CFLAGS)

pointingangles.o: pointingangles.c
	$(CC) -c pointingangles.c -o pointingangles.o $(CFLAGS)

plotcircles.o: plotcircles.c
	$(CC) -c plotcircles.c -o plotcircles.o $(CFLAGS)

plotting_copol.o: plotting_copol.c
	$(CC) -c plotting_copol.c -o plotting_copol.o $(CFLAGS)

plotting_crosspol.o: plotting_crosspol.c
	$(CC) -c plotting_crosspol.c -o plotting_crosspol.o $(CFLAGS)

main.o: main.c
	$(CC) -c main.c -o main.o $(CFLAGS)

getarrays.o: getarrays.c
	$(CC) -c getarrays.c -o getarrays.o $(CFLAGS)

fitphase.o: fitphase.c
	$(CC) -c fitphase.c -o fitphase.o $(CFLAGS)

efficiency.o: efficiency.c
	$(CC) -c efficiency.c -o efficiency.o $(CFLAGS)

outputfilefunctions.o: outputfilefunctions.c
	$(CC) -c outputfilefunctions.c -o outputfilefunctions.o $(CFLAGS)

SCANDATA.o: SCANDATA.c
	$(CC) -c SCANDATA.c -o SCANDATA.o $(CFLAGS)

utilities.o: utilities.c
	$(CC) -c utilities.c -o utilities.o $(CFLAGS)

z.o: z.c
	$(CC) -c z.c -o z.o $(CFLAGS)
