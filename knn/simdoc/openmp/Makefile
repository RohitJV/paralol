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
	./sd -nnbrs=50 -minsim=.01 -nthread=1 /export/scratch/CSCI5451_S18/assignment-3/sports.mat > out50-1
	./sd -nnbrs=50 -minsim=.01 -nthread=2 /export/scratch/CSCI5451_S18/assignment-3/sports.mat > out50-2
	./sd -nnbrs=50 -minsim=.01 -nthread=4 /export/scratch/CSCI5451_S18/assignment-3/sports.mat > out50-4
	./sd -nnbrs=50 -minsim=.01 -nthread=8 /export/scratch/CSCI5451_S18/assignment-3/sports.mat > out50-8
	./sd -nnbrs=50 -minsim=.01 -nthread=1 /export/scratch/CSCI5451_S18/assignment-3/sports5.mat > out50-5-1
	./sd -nnbrs=50 -minsim=.01 -nthread=2 /export/scratch/CSCI5451_S18/assignment-3/sports5.mat > out50-5-2
	./sd -nnbrs=50 -minsim=.01 -nthread=4 /export/scratch/CSCI5451_S18/assignment-3/sports5.mat > out50-5-4
	./sd -nnbrs=50 -minsim=.01 -nthread=8 /export/scratch/CSCI5451_S18/assignment-3/sports5.mat > out50-5-8
	./sd -nnbrs=50 -minsim=.01 -nthread=1 /export/scratch/CSCI5451_S18/assignment-3/sports10.mat > out50-10-1
	./sd -nnbrs=50 -minsim=.01 -nthread=2 /export/scratch/CSCI5451_S18/assignment-3/sports10.mat > out50-10-2
	./sd -nnbrs=50 -minsim=.01 -nthread=4 /export/scratch/CSCI5451_S18/assignment-3/sports10.mat > out50-10-4
	./sd -nnbrs=50 -minsim=.01 -nthread=8 /export/scratch/CSCI5451_S18/assignment-3/sports10.mat > out50-10-8
	./sd -nnbrs=50 -minsim=.01 -nthread=1 /export/scratch/CSCI5451_S18/assignment-3/sports20.mat > out50-20-1
	./sd -nnbrs=50 -minsim=.01 -nthread=2 /export/scratch/CSCI5451_S18/assignment-3/sports20.mat > out50-20-2
	./sd -nnbrs=50 -minsim=.01 -nthread=4 /export/scratch/CSCI5451_S18/assignment-3/sports20.mat > out50-20-4
	./sd -nnbrs=50 -minsim=.01 -nthread=8 /export/scratch/CSCI5451_S18/assignment-3/sports20.mat > out50-20-8
	./sd -nnbrs=50 -minsim=.01 -nthread=1 /export/scratch/CSCI5451_S18/assignment-3/sports30.mat > out50-30-1
	./sd -nnbrs=50 -minsim=.01 -nthread=2 /export/scratch/CSCI5451_S18/assignment-3/sports30.mat > out50-30-2
	./sd -nnbrs=50 -minsim=.01 -nthread=4 /export/scratch/CSCI5451_S18/assignment-3/sports30.mat > out50-30-4
	./sd -nnbrs=50 -minsim=.01 -nthread=8 /export/scratch/CSCI5451_S18/assignment-3/sports30.mat > out50-30-8 
	./sd -nnbrs=200 -minsim=.01 -nthread=1 /export/scratch/CSCI5451_S18/assignment-3/sports.mat > out200-1
	./sd -nnbrs=200 -minsim=.01 -nthread=2 /export/scratch/CSCI5451_S18/assignment-3/sports.mat > out200-2
	./sd -nnbrs=200 -minsim=.01 -nthread=4 /export/scratch/CSCI5451_S18/assignment-3/sports.mat > out200-4
	./sd -nnbrs=200 -minsim=.01 -nthread=8 /export/scratch/CSCI5451_S18/assignment-3/sports.mat > out200-8
	./sd -nnbrs=200 -minsim=.01 -nthread=1 /export/scratch/CSCI5451_S18/assignment-3/sports5.mat > out200-5-1
	./sd -nnbrs=200 -minsim=.01 -nthread=2 /export/scratch/CSCI5451_S18/assignment-3/sports5.mat > out200-5-2
	./sd -nnbrs=200 -minsim=.01 -nthread=4 /export/scratch/CSCI5451_S18/assignment-3/sports5.mat > out200-5-4
	./sd -nnbrs=200 -minsim=.01 -nthread=8 /export/scratch/CSCI5451_S18/assignment-3/sports5.mat > out200-5-8
	./sd -nnbrs=200 -minsim=.01 -nthread=1 /export/scratch/CSCI5451_S18/assignment-3/sports10.mat > out200-10-1
	./sd -nnbrs=200 -minsim=.01 -nthread=2 /export/scratch/CSCI5451_S18/assignment-3/sports10.mat > out200-10-2
	./sd -nnbrs=200 -minsim=.01 -nthread=4 /export/scratch/CSCI5451_S18/assignment-3/sports10.mat > out200-10-4
	./sd -nnbrs=200 -minsim=.01 -nthread=8 /export/scratch/CSCI5451_S18/assignment-3/sports10.mat > out200-10-8
	./sd -nnbrs=200 -minsim=.01 -nthread=1 /export/scratch/CSCI5451_S18/assignment-3/sports20.mat > out200-20-1
	./sd -nnbrs=200 -minsim=.01 -nthread=2 /export/scratch/CSCI5451_S18/assignment-3/sports20.mat > out200-20-2
	./sd -nnbrs=200 -minsim=.01 -nthread=4 /export/scratch/CSCI5451_S18/assignment-3/sports20.mat > out200-20-4
	./sd -nnbrs=200 -minsim=.01 -nthread=8 /export/scratch/CSCI5451_S18/assignment-3/sports20.mat > out200-20-8
	./sd -nnbrs=200 -minsim=.01 -nthread=1 /export/scratch/CSCI5451_S18/assignment-3/sports30.mat > out200-30-1
	./sd -nnbrs=200 -minsim=.01 -nthread=2 /export/scratch/CSCI5451_S18/assignment-3/sports30.mat > out200-30-2
	./sd -nnbrs=200 -minsim=.01 -nthread=4 /export/scratch/CSCI5451_S18/assignment-3/sports30.mat > out200-30-4
	./sd -nnbrs=200 -minsim=.01 -nthread=8 /export/scratch/CSCI5451_S18/assignment-3/sports30.mat > out200-30-8 