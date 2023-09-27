#include "point.h"

#ifndef CELL_H
#define CELL_H

extern const double CELL_SIDE_LENGTH;

typedef struct
{
    Point bottom_left;
    Point bottom_right;
    Point top_left;
    Point top_right;
} Cell;

void create_cell(Cell *cell, double bottom_left_x, double bottom_left_y);

void print_cell(Cell *cell);

#endif