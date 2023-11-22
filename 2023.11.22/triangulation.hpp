#pragma once
#include "header.hpp"
#include "infrastructure.hpp"

// double **gift_wrapping_2d(int num_segments, int num_seeds, double error);
double **gift_wrapping_2d(double **segments, double **vertices, double error);

double *positively_oriented_triangle_2d(double p1_x, double p1_y, double p2_x, double p2_y, double **vertices, double** segments, double error);

double orientation_2d(double a_x, double a_y, double b_x, double b_y, double c_x, double c_y);

double in_circle(double a_x, double a_y, double b_x, double b_y, double c_x, double c_y, double d_x, double d_y);

double triple_product(double a_x, double a_y, double a_z, double b_x, double b_y, double b_z, double c_x, double c_y, double c_z);

double if_itersection(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y, double p4_x, double p4_y, double error);