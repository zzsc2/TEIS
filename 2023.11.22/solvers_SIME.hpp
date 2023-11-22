#include "header.hpp"
#include "integrals.hpp"
#include "infrastructure.hpp"

double *FuncJperpendicular(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***B_tri_3_18, double ***pressure_2_tri_18, int position_sec, int position_tri);
double *FuncGradSigma(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****Jperpendicular_2_tri_3_18, double ***B_tri_3_18, int position_sec, int position_tri);
double FuncDivGradSigma(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****GradSigma_2_tri_3_18, int position_sec, int position_tri);
double *FuncJparallel(double u, double v, double ***G, double ***B_tri_3_18, double **sigma_tri_18, int position_sec, int position_tri);
double *FuncCurlJ(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****Jperpendicular_2_tri_3_18, int position_sec, int position_tri);



double *FuncGradP(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***pressure_2_tri_18, int position_sec, int position_tri);
double *FuncGradPparallel(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***B_tri_3_18, double ***gradP_tri_3_18, int position_sec, int position_tri);
double FuncDivGradPparallel(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****GradPparallel_2_tri_3_18, int position_sec, int position_tri);



double FuncPressure(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, double ***pressure_3_tri_18, double ****GradPparallel_2_tri_3_18, double ***B_tri_3_18, int position_sec, int position_tri);


double ***solve_magnetic_field_sparse_LM(double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangels_in_plane, int num_scalar, double error, time_t time_s);
double ***solve_pressure_sparse_LM(double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangels_in_plane, double error, time_t time_s);

double pressure_initialization(double u, double v, int position_triangle, double **triangles_a_b_c_theta_x0_y0, vector<double> coefficients);








