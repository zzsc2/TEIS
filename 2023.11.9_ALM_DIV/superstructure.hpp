#pragma once
#include "header.hpp"
#include "solvers.hpp"
#include "infrastructure.hpp"

double *biot_savart(double r, double z, double phi, double ***coils_p);

double ***divergence_cleaning(double ***B_3_tri_18N, double ***Mpoisson_tri_18N_18N, int nfp_ft, int num_modes, int num_sections, \
double ***G, double **GLs, int GL_order, double ***R_sec_tri_18, double ***Z_sec_tri_18, double ***Rzeta_sec_tri_18, double ***Zzeta_sec_tri_18, double **triangles_uv, time_t time_s, \
double ***M1, \
double *****MG, double ****M2_integral4, double ***BG_tri_18_18, double error, \
int num_edges, int **points_sequence, double **boundary_DoFs);



