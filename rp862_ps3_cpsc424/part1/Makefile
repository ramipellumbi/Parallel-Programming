#
# This Makefile assumes that the following module files are loaded:
#
#   iomkl
#
# Set the compiler command
CC = mpicc

# Location of timing routines

TIMING = /gpfs/gibbs/project/cpsc424/shared/utils/timing

# Set the compiler options

CFLAGS = -g -O3 -xHost -fno-alias -std=c99 -I$(TIMING)

all: step1 

task1: task1.o $(TIMING)/timing.o rwork.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

# Be careful not to remove rwork.o. If you do, you'll need to get a new copy.
clean:
	rm -f task1 task1.o

