#include "header.hpp"
#include "infrastructure.hpp"
#include "triangulation.hpp"

double **ruppert_refinement(int num_segments, int num_seeds, double range_of_seeds, double rho, double error);

double **chew_first_refinement(int num_segments, double rho, double scale, double error);

double **split_segment(double *segment, double **segments);

double *find_circumcenter(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y);

double in_min_containment_circle(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y, double error);

double radius_edge_ratio(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y);

double length_of_radius(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y);

double length_of_min_edge(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y);

double *adjacent(double x1, double y1, double x2, double y2, double** triangles, double **vertices, double **segments, double error);

double **dig_cavity(double x1, double y1, double x2, double y2, double x3, double y3, double** triangles, double **vertices, double **segments, double error);

double **bowyer_watson(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double** triangles, double **vertices, double **segments, double error);