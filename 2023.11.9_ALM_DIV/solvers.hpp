#include "header.hpp"
#include "infrastructure.hpp"

double *FuncPoisson(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri);
// test
double ***FuncBeltrami(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri, int num_scalar, double mu);
double **boundary_DoFs_beltrami(double ***triangles_RZ, int **points_sequence, double ***coils_p, int num_sections, int num_modes, int nfp_ft);
double FuncDivergence(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ****B_2_tri_3_18);
double *FuncProjection(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ***B_tri_3_18, double ***poisoon_2_tri_18);
double ***FuncProjection2(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ***B_tri_3_18, double ***poisoon_2_tri_18);
double *FuncDivCoeff(int position_p, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, int position_sec, int position_tri, int num_sections);

double ****divergence_free_constraints(double **triangles_uv, int num_sections, int num_modes, int nfp_ft, double ***R, double ***Z, double ***Rzeta, double ***Zzeta);
double *****MG_tri_6_18_18_18(double ***G, double **triangles_uv);
double ***BG_tri_18_18(double ***G, double **triangles_uv);
double ****M_integral4(int nfp, int num_modes, double error);
double ***M_tri_18Nn_18Nn(double *****MG_tri_6_18_18_18, double *****object_tri_18N_ope_3_3, double ****array_integral4, int num_triangles_in_plane, int num_modes, int nfp, int num_scalar, double error, int num_derivates);
double ***Mpoisson_tri_18N_18N(double *****MG_tri_6_18_18_18, double ***poisson_tri_18N_10, double ****array_integral4, int num_triangles_in_plane, int num_modes, int nfp, int num_scalar, double error);
double **B_tri_18Nn(double ***BG_tri_18_18, double *****object_tri_18N_ope_3_3, int num_triangles_in_plane, int num_modes, int nfp, int num_scalar, double error, int num_derivates);
double ***solve_magnetic_field_sparse_direct(double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangels_in_plane, int num_scalar, double error, time_t time_s);
double ***solve_magnetic_field_sparse_guess(double ***guess_sca_tri_18N, double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangels_in_plane, int num_scalar, double error, time_t time_s);
double ***solve_magnetic_field_sparse_lagrange_multiplier(double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangels_in_plane, int num_scalar, int flag, double error, time_t time_s);
double ***solve_magnetic_field_sparse_lagrange_multiplier_div(double**** divergence_tri_3_var_modes, double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangels_in_plane, int nfp_ft, int num_scalar, int flag, double error, time_t time_s);
double ***solve_magnetic_field_sparse_ALM_div(double**** divergence_tri_3_var_modes, double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangels_in_plane, int nfp_ft, int num_scalar, int flag, double error, time_t time_s);

double pressure_initialization(double u, double v, int position_triangle, double **triangles_a_b_c_theta_x0_y0, vector<double> coefficients);
double ***Func_SIME_magnetic_field(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***mf_3_tri_18, double **pres_tri_18, double ***G, int num_triangles, int num_modes, int nfp_ft, int position_sec, int position_tri, int num_scalar);





