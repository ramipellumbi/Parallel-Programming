#include "point.h"

const double CELL_SIDE_LENGTH = 0.001;

/*
 * Add two points as complex number
 */
Point add_points(Point p1, Point p2)
{
    Point sum;
    sum.x = p1.x + p2.x;
    sum.y = p1.y + p2.y;

    return sum;
}

double magnitude_of_point_squared(Point point)
{
    return point.x * point.x + point.y * point.y;
}

/**
 * Multiply two points as complex numbers
 */
Point multiply_points(Point p1, Point p2)
{
    Point product;
    product.x = p1.x * p2.x - p1.y * p2.y; // Real part
    product.y = p1.x * p2.y + p1.y * p2.x; // Imaginary part
    return product;
}

void print_point(Point p)
{
    printf("Point: (x = %.2f, y = %.2f)\n", p.x, p.y);
}