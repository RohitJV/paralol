CC = gcc 
OPTFLAGS = -O3 -fopenmp -fno-strict-aliasing -D_GNU_SOURCE
COPTIONS = -DLINUX -D_FILE_OFFSET_BITS=64 -std=c99 -Wall\
           -Wno-unused-function -Wno-unused-label -Wno-unused-variable\
           -Wno-parentheses -Wsequence-point

#Standard Libraries
STDLIBS = -lm

CFLAGS = $(COPTIONS)  $(OPTFLAGS)
LIBS = $(STDLIBS)

SOURCE = sd.c getopt.c cmdline.c util.c
EXE = sd

default:
	export OMP_PLACES=cores
	export OMP_PROC_BIND=spread
	$(CC) $(CFLAGS) $(SOURCE) $(LIBS) -o $(EXE)
