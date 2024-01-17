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

/**
 * Leapfrog Method
 *
 * 1. Each thread needs to be initialized with a seed that is a function of the master seed and the thread's ID.
 * 2. The thread updates its seed by advancing the seed by the number of threads.
 *
 * Thread 0 uses numbers at position 0, NUM_THREADS, 2*NUM_THREADS,...
 * Thread 1 uses numbers at position 1, NUM_THREADS + 1, 2*NUM_THREADS + 1,...
 * ...
 */

// Advance the RNG state by `n` steps
void leapfrog(uint64_t *seed, uint64_t n)
{
    for (uint64_t i = 0; i < n; ++i)
    {
        *seed = 6364136223846793005ULL * (*seed) + 1;
    }
}

void dsrand_parallel_ts(unsigned s)
{
    int thread_id = omp_get_thread_num();

    uint64_t master_seed = s - 1;
    seed_ts = master_seed;

    // Each thread leapfrogs `thread_id` steps from master_seed
    leapfrog(&seed_ts, thread_id);
}

double drand_parallel_ts()
{
    // leapfrog the seed ahead num_threads
    int num_threads = omp_get_num_threads();
    leapfrog(&seed_ts, num_threads);

    return ((double)(seed_ts >> 33) / (double)RAND_MAX);
}
