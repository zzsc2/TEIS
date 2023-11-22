#pragma once
#include "header.hpp"
#include "infrastructure.hpp"

int m_order[20];
int n_order[20];
void mn_order();

double integral1(int M, int N, double a, double b, double c);

double integral2(double a, double b, double c, double ***g, int i, int j, int num);

double integral3(double ***G, int i, int j, int k, int position, int m, int n, double *****array_integral3);

double integral4(int nfp, int ntor, int t, int k, int j, int Gj, double error);

double **Gauss_Legendre_quadrature_5();

double numerical_integration_over_triangle(double A, double B, double C, std::function<double(double, double)> func, double **GLs, int n); 
// test
double numerical_integration1_over_triangle(double A, double B, double C, std::function<double(double, double, double, double, double, double, double**, int, double, double, int, int)> func, \
double theta, double x0, double y0, double zeta, double **fourier_coefficients, int nfp, double sine, double cosine, int position, int column, \
double **GLs, int n);
// debug


//
double **numerical_integration4_over_triangle(double A, double B, double C, std::function<double**(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int, int, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri, int column, int rank, \
double **GLs, int n);


double ***integral_functional1_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double(double, double, double, double, double, double, double **, int, double, double)> functional, \
int nfp, double **fourier_coefficients, double **triangles_uv);
// debug

double *****integral_functional4_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double **(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double **triangles_uv, int rank);
// P

// divergence-free constraint
double *numerical_integration_DivCoeff_over_triangle(double A, double B, double C, std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int column, \
double **GLs, int n);

double ****integral_DivCoeff_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double **triangles_uv);

// projection
double *numerical_integration_projection_over_triangle(double A, double B, double C, std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int, int, int, double ***, double ***, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ***B_3_tri_18N, double ***poisoon_1_tri_18N, int column, \
double **GLs, int n);

double ****integral_projection_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int, int, int, double ***, double ***)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double **triangles_uv, int num_modes, int nfp_ft, double ***B_3_tri_18N, double ***poisoon_1_tri_18N);

// divergence
double numerical_integration_divergence_over_triangle(double A, double B, double C, std::function<double(double , double , double ***, double ***, double ***, double ***, double ***, int , int , int , int , int, double ****, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ****B_2_tri_3_18, int column, \
double **GLs, int n);

double ***integral_divergence_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double (double , double , double ***, double ***, double ***, double ***, double ***, int , int , int , int , int, double ****)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double **triangles_uv, int num_modes, int nfp_ft, double ***B_3_tri_18N);

// poisson
double *numerical_integration_poisson_over_triangle(double A, double B, double C, std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int, int, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri, int column, int num_variables, \
double **GLs, int n);

double ****integral_poisson_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double **triangles_uv, int num_variables);
// 

// R Z Rzeta Zzeta Rzetazeta Zzetazeta
double *numerical_integration_RZ_over_triangle(double A, double B, double C, std::function<double *(double, double, double, double, double, double**, double**, double**, int, int, double, double, int, int)> func, \
double x0, double y0, double zeta, double **fc, double **fc_zeta, double **fc_zetazeta, int num_variables, int nfp, double sine, double cosine, int position, int column, \
double **GLs, int n);

double ****integral_RZ_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double *(double, double, double, double, double, double **, double **, double **, int, int, double, double)> functional, \
int nfp, double **fc, double **fc_zeta, double **fc_zetazeta, int num_variables, double **triangles_uv);
// 

// jacobian
double numerical_integration_jacobian_over_triangle(double A, double B, double C, std::function<double(double, double, double ***, double ***, double ***, int, int, int)> func, \
double ***R, double ***Z, double ***G, int position_sec, int position_tri, int column, \
double **GLs, int n);

double ***integral_jacobian_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double(double u, double v, double ***R, double ***Z, double ***G, int position_sec, int position_tri)> functional, \
double ***R, double ***Z, double **triangles_uv);
// 

// beltrami
double ***numerical_integration_beltrami_over_triangle(double A, double B, double C, std::function<double***(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int, int, double, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri, int num_scalar, double mu, \
int column, \
double **GLs, int n);

double ******integral_beltrami_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double ***(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int, int, double)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double **triangles_uv, int num_scalar, double mu);
//

// SMIE_magnetic_field
double ***numerical_integration_SMIEmf_over_triangle(double A, double B, double C, std::function<double***(double , double , double ***, double ***, double ***, double ***, double ***, double **, double ***, int , int , int , int , int , int, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***mf_3_tri_18, double **pres_tri_18, double ***G, int num_triangles, int num_modes, int nfp_ft, int position_sec, int position_tri, int num_scalar, \
int column, \
double **GLs, int n);

double ******integral_SMIEmf_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double ***(double , double , double ***, double ***, double ***, double ***, double ***, double **, double ***, int , int , int , int , int , int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***mf_3_tri_18N, double **pres_tri_18N, double **triangles_uv, int num_modes, int nfp_ft, int num_scalar);
//






