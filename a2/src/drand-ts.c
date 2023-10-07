#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <omp.h>

// using thread_local from <threads.h> instead of static had
// same impact as threadprivate
static uint64_t seed_ts;

// add for parallel fix - tracks if a thread has been assigned a seed
static bool is_thread_seed_initialized = false;

// ensures each thread is modifying it's own seed
#pragma omp threadprivate(seed_ts, is_thread_seed_initialized)

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

void dsrand_parallel_ts(unsigned s)
{
    seed_ts = s - 1;
    is_thread_seed_initialized = true;
    printf("Seed_ts = %lu. RAND_MAX = %d.\n", seed_ts, RAND_MAX);
}

double drand_parallel_ts()
{
    // if the thread has not initialized it's own unique seed do so
    if (!is_thread_seed_initialized)
    {
        int offset = 10;
        int prime = 23;
        int thread_id = omp_get_thread_num();
        // ensure this thread gets a unique seed
        uint64_t thread_seed = offset + thread_id * prime;
        dsrand_ts(thread_seed);
    }

    seed_ts = 6364136223846793005ULL * seed_ts + 1;
    return ((double)(seed_ts >> 33) / (double)RAND_MAX);
}