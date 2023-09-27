#include "cell.h"
#include "point.h"
#include <stdio.h>

void create_cell(Cell *cell, double bottom_left_x, double bottom_left_y)
{
    // Update bottom_left point
    cell->bottom_left.x = bottom_left_x;
    cell->bottom_left.y = bottom_left_y;

    // Update bottom_right point
    cell->bottom_right.x = bottom_left_x + CELL_SIDE_LENGTH;
    cell->bottom_right.y = bottom_left_y;

    // Update top_left point
    cell->top_left.x = bottom_left_x;
    cell->top_left.y = bottom_left_y + CELL_SIDE_LENGTH;

    // Update top_right point
    cell->top_right.x = bottom_left_x + CELL_SIDE_LENGTH;
    cell->top_right.y = bottom_left_y + CELL_SIDE_LENGTH;
}

void print_cell(Cell *cell)
{
    printf("Cell:\n");
    printf("  Bottom Left: ");
    print_point(cell->bottom_left);
    printf("  Bottom Right: ");
    print_point(cell->bottom_right);
    printf("  Top Left: ");
    print_point(cell->top_left);
    printf("  Top Right: ");
    print_point(cell->top_right);
}