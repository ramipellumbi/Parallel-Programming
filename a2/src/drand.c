#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <threads.h>

#include "drand.h"

static uint64_t seed; 

void dsrand(unsigned s)
{
    seed = s - 1;
    printf("Seed = %lu. RAND_MAX = %d.\n", seed, RAND_MAX);
}

double drand(void)
{
    seed = 6364136223846793005ULL * seed + 1;
    return ((double)(seed >> 33) / (double)RAND_MAX);
}

// using thread_local much faster than #omp critical section
// https://stackoverflow.com/questions/71609066/are-thread-local-objects-initialized-to-0-in-c
// https://gcc.gnu.org/onlinedocs/gcc/Thread-Local.html
thread_local uint64_t seed_ts;

void dsrand_ts(unsigned s)
{
    seed_ts = s - 1;
    printf("Seed_ts = %lu. RAND_MAX = %d.\n", seed_ts, RAND_MAX);
}

double drand_ts(void)
{
    seed_ts = 6364136223846793005ULL * seed_ts + 1;
    return ((double)(seed_ts >> 33) / (double)RAND_MAX);
}