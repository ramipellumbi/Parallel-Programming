#include "drand.h"

#ifndef UTILITIES_H
#define UTILITIES_H

static inline double get_random_double_in_bounds(double min, double max)
{
    return min + drand() * (max - min);
}

static inline double get_random_double_in_bounds_ts(double min, double max)
{
    return min + drand_ts() * (max - min);
}

int get_environment_value(const char *env_name);

#endif