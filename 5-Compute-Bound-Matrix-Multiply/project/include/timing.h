#ifndef TIMING_H
#define TIMING_H

#include <sys/time.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/resource.h>

void timing(double *wcTime, double *cpuTime);

#endif