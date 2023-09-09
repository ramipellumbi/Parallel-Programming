# Assumes that you have already loaded the intel compiler:
#        module load intel
#
# Assumes that you've copied timing.c and timing.h to this directory

CC = icc

#Various compiler options

#CFLAGS = -g -O0 -fno-alias -std=c99
#CFLAGS = -g -O1 -fno-alias -std=c99
#CFLAGS = -g -O0 -no-vec -no-simd -fno-alias -std=c99 -fp-model precise
#CFLAGS = -O0 -no-vec -no-simd -fno-alias -std=c99
#CFLAGS = -O1 -no-vec -no-simd -fno-alias -std=c99
#CFLAGS = -O3 -no-vec -no-simd -fno-alias -std=c99
#CFLAGS = -g -O3 -xHost -fno-alias -std=c99 -qopt-report=5 -fsource-asm
#CFLAGS = -std=c99

CFLAGS = -g -O3 -xHost -fno-alias -std=c99

# Build code for Exercise 1
pi:	pi.o timing.o
	$(CC) -o $@ $(CFLAGS) $^

# Build code for Exercise 2
# Note: dummy.c is provided to you
triad:  triad.o timing.o dummy.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f pi triad *.o
