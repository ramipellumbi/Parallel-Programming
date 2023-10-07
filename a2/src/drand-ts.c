#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <omp.h>

// using thread_local from <threads.h> instead of static had
// same impact as threadprivate
static uint64_t seed_ts;

// ensures each thread is modifying it's own seed
#pragma omp threadprivate(seed_ts)

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