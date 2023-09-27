#include "drand.h"
#include "cell.h"
#include "point.h"

double get_random_double_in_bounds(double min, double max)
{
    return min + drand() * (max - min);
}

Point generate_random_point_in_cell(Cell *cell)
{
    double x = get_random_double_in_bounds(cell->bottom_left.x, cell->bottom_right.x);
    double y = get_random_double_in_bounds(cell->bottom_left.y, cell->top_left.y);

    return (Point){.x = x, .y = y};
}