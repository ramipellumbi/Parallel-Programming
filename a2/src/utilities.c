#include <stdio.h>
#include <stdlib.h>
#include "utilities.h"

int get_environment_value(const char *env_name)
{
    char *value_str = getenv(env_name);
    if (value_str == NULL)
    {
        return -1;
    }
    int value = atoi(value_str);

    return value;
}