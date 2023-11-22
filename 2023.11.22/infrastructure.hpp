#pragma once
#include "header.hpp"

// #define pi 3.14159265358979
// #define mu0 0.000001256637061435917
// #define bs mu0/(4*pi)
// #define NT 5

void timing(time_t time_s);
void timing();

template<class T>
T power_nonnegative(T x, int n);

double factorial(double n);

int n_table = 100000;
double step_r = 2*pi / n_table;
double *cosine_table = new double [n_table];
double *sine_table = new double [n_table];
void table();

double cosine_lookup_approximation(double angle);

double sine_lookup_approximation(double angle);

double cosine_MMP_approximation(double angle);

double sine_MMP_approximation(double angle);


int binomial_coefficient(int n, int k);

template<class T>
bool value_comparison(T a, T b, double error);

template<class T>
T *enlarge_1d(T *old_array, int old_size, int new_size);

template<class T>
T *reduce_1d(T *old_array, int old_size, int location, int reduced_size);

template<class T>
T **reduce_2d(T **old_array, int old_size, int location, int reduced_size, int dimension);

template<class T>
T **enlarge_2d(T **old_array, int old_size, int new_size, int dimension);

template<class T>
T **insert_2d(T **old_array, int old_size, int location, int added_size, int dimension);

double fourier_basis(int i, int nfp, double zeta);
double fourier_basis_test(int i, int nfp, double zeta);
double fourier_basis_test2(int i, int nfp, double zeta);
double fourier_basis_test3(int i, int nfp, double zeta);

void print_RZ(double **R, double **Z, double ***G, double **triangles_uv, int num_modes, int nfp_ft);
void print_B(double ***F, double ***G, double **triangles_uv, int num_modes, int nfp_ft, int num_scalar);
void print_p(double **F, double ***G, double **triangles_uv, int num_modes, int nfp_ft);

double low_of_cosines(double a, double b, double c);
