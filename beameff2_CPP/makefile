# Project: beameff2_CPP
# Makefile for Linux only
# Build windows executable with Eclipse internal builder for MinGW/GCC 3.4.5

CPP  = g++ -DLINUX -g
CC   = gcc -DLINUX -g
OBJ  = iniparser4.0/dictionary.o iniparser4.0/iniparser.o unwrap2D/unwrap2D.o \
	nr/brent.o nr/f1dim.o nr/mnbrak.o nr/linmin.o nr/frprmn.o nr/nrutil.o \
	main.o ALMAConstants.o AnalyzeResults.o BeamEffInputFile.o EfficiencyData.o \
	FitAmplitude.o FitPhase.o ScanDataRaster.o ScanData.o ScanSet.o SWVersion.o
LIBS = -lm
INC = -I./iniparser4.0 -I./nr -I./util -I./unwrap2D
BIN  = beameff2_64
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before beameff2_64 all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(OBJ) -o "beameff2_64" $(LIBS)

dictionary.o: iniparser4.0/dictionary.c
	$(CC) -c iniparser4.0/dictionary.c -o iniparser4.0/dictionary.o $(INC)

iniparser.o: iniparser4.0/iniparser.c
	$(CC) -c iniparser4.0/iniparser.c -o iniparser4.0/iniparser.o $(INC)

unwrap2D.o: unwrap2D/unwrap2D.c
	$(CC) -c unwrap2D/unwrap2D.c -o unwrap2D/unwrap2D.o $(INC)

brent.o: nr/brent.c
	$(CC) -c nr/brent.c -o nr/brent.o $(INC)

f1dim.o: nr/f1dim.c
	$(CC) -c nr/f1dim.c -o nr/f1dim.o $(INC)

mnbrak.o: nr/mnbrak.c
	$(CC) -c nr/mnbrak.c -o nr/mnbrak.o $(INC)

linmin.o: nr/linmin.c
	$(CC) -c nr/linmin.c -o nr/linmin.o $(INC)

frprmn.o: nr/frprmn.c
	$(CC) -c nr/frprmn.c -o nr/frprmn.o $(INC)

nrutil.o: nr/nrutil.c
	$(CC) -c nr/nrutil.c -o nr/nrutil.o $(INC)

main.o: main.cpp
	$(CPP) -c main.cpp -o main.o $(INC)

ALMAConstants.o: ALMAConstants.cpp
	$(CPP) -c ALMAConstants.cpp -o ALMAConstants.o $(INC)
	
AnalyzeResults.o: AnalyzeResults.cpp
	$(CPP) -c AnalyzeResults.cpp -o AnalyzeResults.o $(INC)
	
BeamEffInputFile.o: BeamEffInputFile.cpp
	$(CPP) -c BeamEffInputFile.cpp -o BeamEffInputFile.o $(INC)

EfficiencyData.o: EfficiencyData.cpp
	$(CPP) -c EfficiencyData.cpp -o EfficiencyData.o $(INC)

FitAmplitude.o: FitAmplitude.cpp
	$(CPP) -c FitAmplitude.cpp -o FitAmplitude.o $(INC)

FitPhase.o: FitPhase.cpp
	$(CPP) -c FitPhase.cpp -o FitPhase.o $(INC)

ScanDataRaster.o: ScanDataRaster.cpp
	$(CPP) -c ScanDataRaster.cpp -o ScanDataRaster.o $(INC)

ScanData.o: ScanData.cpp
	$(CPP) -c ScanData.cpp -o ScanData.o $(INC)

ScanSet.o: ScanSet.cpp
	$(CPP) -c ScanSet.cpp -o ScanSet.o $(INC)

SWVersion.o: SWVersion.cpp
	$(CPP) -c SWVersion.cpp -o SWVersion.o $(INC)
