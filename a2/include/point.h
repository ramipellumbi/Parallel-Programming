#ifndef POINT_H
#define POINT_H

typedef struct
{
    double x;
    double y;
} Point;

Point add_points(Point p1, Point p2);

double magnitude_of_point_squared(Point point);

Point multiply_points(Point p1, Point p2);

void print_point(Point p);

#endif