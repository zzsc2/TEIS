#include "header.hpp"
#include "infrastructure.hpp"



// pressure
double numerical_integration_pressure_over_triangle(double A, double B, double C, std::function<double(double , double , int , double **, vector<double>, int, double ***)> func, \
int position_tri, double **triangles_uv, vector<double> coefficients, int column, double ***G, \
double **GLs, int n);

double ***integral_pressure_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double(double u, double v, int position_tr, double **triangles_uv, vector<double> coefficients)> functional, \
double **triangles_uv, vector<double> coefficients);
//

// projection2
double ***numerical_integration_projection2_over_triangle(double A, double B, double C, std::function<double***(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int, int, int, double ***, double ***, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ***B_tri_3_18, double ***poisoon_2_tri_18, int column, \
double **GLs, int n);

double ******integral_projection2_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double ***(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int, int, int, double ***, double ***)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double **triangles_uv, int num_modes, int nfp_ft, double ***B_3_tri_18N, double ***poisoon_1_tri_18N);

// SIME grad(P)
double *numerical_integration_gradP_over_triangle(double A, double B, double C, std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ***, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***pres_2_tri_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n);

double ****integral_gradP_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double **pres_tri_18N, double **triangles_uv);

// SIME grad(P) parallel
double *numerical_integration_gradPparallel_over_triangle(double A, double B, double C, std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ***, double ***, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***mf_tri_3_18, double ***gradP_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n);

double ****integral_gradPparallel_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***mf_3_tri_18N, double ***gradP_3_tri_18N, double **triangles_uv) ;

// SIME J perpendicular
double *numerical_integration_Jperpendicular_over_triangle(double A, double B, double C, std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ***, double ***, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***mf_tri_3_18, double ***gradP_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n);

double ****integral_Jperpendicular_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***mf_3_tri_18N, double **gradP_3_tri_18N, double **triangles_uv);

// SIME curl(J)
double *numerical_integration_curlJ_over_triangle(double A, double B, double C, std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ****, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****J_2_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n);

double ****integral_curlJ_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ****, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***J_3_tri_18N, double **triangles_uv);

// SIME grad(sigma)
double *numerical_integration_gradSigma_over_triangle(double A, double B, double C, std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, double ****, double ***, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****J_2_tri_3_18, double ***B_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n);

double ****integral_gradSigma_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ****, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***J_3_tri_18N, double ***B_3_tri_18N, double **triangles_uv);

// SIME div(grad(sigma))
double numerical_integration_divGradSigma_over_triangle(double A, double B, double C, std::function<double(double, double, double ***, double ***, double ***, double ***, double ***, double ****, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****GradSigma_2_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n);

double ***integral_divGradSigma_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double (double, double, double ***, double ***, double ***, double ***, double ***, double ****, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***GradSigma_3_tri_18N, double **triangles_uv);

// SIME J parallel
double *numerical_integration_Jparallel_over_triangle(double A, double B, double C, std::function<double*(double, double, double ***, double ***, double **, int , int , int)> func, \
double ***G, double ***mf_3_tri_18, double **sigma_tri_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n);

double ****integral_Jparallel_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double **, int, int)> functional, \
double ***mf_3_tri_18N, double **sigma_tri_18N, double **triangles_uv);







// SIME pressure
double numerical_integration_pres_over_triangle(double A, double B, double C, std::function<double(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, double ***, double ****, double ***, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, double ***pressure_3_tri_18, double ****GradPparallel_2_tri_3_18, double ***B_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n);

double ***integral_pres_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double (double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, double ***, double ****, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double **pressure_tri_18N, double ***GradPparallel_tri_3_18N, double ***B_3_tri_18N, double **triangles_uv);









// SIME div(JxB)
double numerical_integration_divJ_over_triangle(double A, double B, double C, std::function<double(double, double, double ***, double ***, double ***, double ***, double ***, double ****, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****J_2_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n);

double ***integral_divJ_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double (double, double, double ***, double ***, double ***, double ***, double ***, double ****, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***J_3_tri_18N, double **triangles_uv);




