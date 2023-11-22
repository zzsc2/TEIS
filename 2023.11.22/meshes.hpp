#include "header.hpp"
#include "infrastructure.hpp"

double ***triangles_mapping(double **triangles_in_circle, double **fourier_coefficients, int num_sections, int nfp);
double coordinate_R(double u, double v, double x0, double y0, double zeta, double **fourier_coefficients, int nfp, double sine, double cosine); 
double coordinate_Z(double u, double v, double x0, double y0, double zeta, double **fourier_coefficients, int nfp, double sine, double cosine); 
double boundary_R(double theta, double zeta_nfp, double **fourier_coefficients);
double boundary_Z(double theta, double zeta_nfp, double **fourier_coefficients);

void triangles_reorder(double ***triangles_RZ, double ** triangles_xy, int num_sections);
int **point_retrieval(double **triangles_xy, int num_edges, int num_points, double error);
double **triangles_a_b_c_theta_x0_y0(double **triangles_xy);
double *foot_of_perpendicular(double x1, double y1, double x2, double y2, double x3, double y3);
double **normal_vector(double **fc_zeta, double **fc_theta, double **triangles_xy, int **points_sequence, int num_sections, int num_modes, int nfp, int nfp_ft);

double ***G_tri_20_18(double **triangles_uv);
double ***M_tri_18_18(double ***G, double **triangles_uv);

double *coordinate_funRZ(double u, double v, double x0, double y0, double zeta, double **fc, double **fc_zeta, double **fc_zetazeta, int num_variables, int nfp, double sine, double cosine);
double jacobian(double u, double v, double ***R, double ***Z, double ***G, int position_sec, int position_tri);
double **P_9_9(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri);

double ***F_sec_tri_18(double ***M_ij, double ***f_i, int num_sections, int num_triangles_in_plane);
// debug
double ****F_sec_tri_18_num(double ***M_ij, double ****f_i, int num_sections, int num_triangles_in_plane, int num_variables);
double *****F_sec_tri_18_m_m(double ***M_ij, double *****f_i, int num_sections, int num_triangles_in_plane, int rank);
double ******F_sec_tri_18_derivates_m_m(double ***M_ij, double ******f_i, int num_sections, int num_triangles_in_plane, int rank, int num_derivates);





double **re_axis(double **triangles_in_circle, double **fourier_coefficients, int nfp, int n, int order);
double ***triangles_transform2(double **triangles_in_circle, double **fourier_coefficients, double **axis_coefficients, int num_sections, int nfp);

double **boundary_recalculation(double **fourier_coefficients, int nfp, int n);
