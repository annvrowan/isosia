# Den compiler der skal bruges
CC=gcc
CCGRENDEL=icc

# remember to add `source /com/bin/intel.sh` to your ~/.bashrc

# Programmet der skal compiles
BIN=../spm

# Kildekoden der skal bruges for at compile
SRC=spm.c

# Brug nedenstaaende til debug/profilering
CFLAGSDBG=-lm -g -pg -fopenmp

# Brug nedenstaaende til "produktion"
CFLAGS=-lm -O3 -fopenmp -ffast-math -funroll-loops
CFLAGSGRENDEL=-O3 -openmp -funroll-loops

# Kan kaldes med `make` eller `make ../spm`
$(BIN): $(SRC)
	$(CC) $(CFLAGS) $< -o $@

# `make debug`
debug: $(BIN)-debug

grendel: $(SRC)
	$(CCGRENDEL) $(CFLAGSGRENDEL) $< -o $(BIN) -lm

$(BIN)-debug: $(SRC)
	$(CC) $(CFLAGSDBG) $< -o $@

# Ryd op med `make clean`
clean:
	$(RM) $(BIN)
	$(RM) *.o

