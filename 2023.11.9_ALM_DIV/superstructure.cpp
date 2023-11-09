#include "superstructure.hpp"

double *biot_savart(double r, double z, double phi, double ***coils_p) {

    int num_coils = coils_p[0][0][4];
    int num_sections = coils_p[0][0][5];

    double *value = new double [3] ();

    double cos_phi = cos(phi), sin_phi = sin(phi);
    double x = r * cos_phi, y = r * sin_phi;
    double sum_x = 0, sum_y = 0, sum_z = 0;
    
    double **array1 = new double *[num_coils], **array2 = new double *[num_coils];
    double **dx = new double *[num_coils], **dy = new double *[num_coils], **dz = new double *[num_coils];
    double **rx = new double *[num_coils], **ry = new double *[num_coils], **rz = new double *[num_coils];
    for (int i = 0; i < num_coils; i++) {
        array1[i] = new double [num_sections-1] ();
        array2[i] = new double [num_sections-1] ();
        dx[i] = new double [num_sections-1] ();
        dy[i] = new double [num_sections-1] ();
        dz[i] = new double [num_sections-1] ();
        rx[i] = new double [num_sections-1] ();
        ry[i] = new double [num_sections-1] ();
        rz[i] = new double [num_sections-1] ();
    }

    #pragma omp parallel for collapse(2) reduction(+:sum_x) reduction(+:sum_y) reduction(+:sum_z)
    for (int i = 0; i < num_coils; i++) {
        for (int j = 0; j < num_sections - 1; j++) {
            dx[i][j] = coils_p[i][j+1][0] - coils_p[i][j][0];
            dy[i][j] = coils_p[i][j+1][1] - coils_p[i][j][1];
            dz[i][j] = coils_p[i][j+1][2] - coils_p[i][j][2];

            rx[i][j] = x - coils_p[i][j][0];
            ry[i][j] = y - coils_p[i][j][1];
            rz[i][j] = z - coils_p[i][j][2];

            array1[i][j] = sqrt(rx[i][j]*rx[i][j] + ry[i][j]*ry[i][j] + rz[i][j]*rz[i][j]);
            array2[i][j] = bs * (coils_p[i][j][3] * 1000000.0) / (array1[i][j] * array1[i][j] * array1[i][j]);

            sum_x += (rz[i][j]*dy[i][j] - ry[i][j]*dz[i][j]) * array2[i][j];
            sum_y += (rx[i][j]*dz[i][j] - rz[i][j]*dx[i][j]) * array2[i][j];
            sum_z += (ry[i][j]*dx[i][j] - rx[i][j]*dy[i][j]) * array2[i][j];
        }
    }

    value[0] =  sum_x * cos_phi + sum_y * sin_phi;
    value[1] =  sum_z;
    value[2] = -sum_x * sin_phi + sum_y * cos_phi;
    
    for (int i = 0; i < num_coils; i++) {
        delete [] dx[i];
        delete [] dy[i];
        delete [] dz[i];
        delete [] rx[i];
        delete [] ry[i];
        delete [] rz[i];
        delete [] array1[i];
        delete [] array2[i];
    }

    delete [] dx;
    delete [] dy;
    delete [] dz;
    delete [] rx;
    delete [] ry;
    delete [] rz;
    delete [] array1;
    delete [] array2;

    return value;
}

double ***divergence_cleaning(double ***B_3_tri_18N, double ***Mpoisson_tri_18N_18N, int nfp_ft, int num_modes, int num_sections, \
double ***G, double **GLs, int GL_order, double ***R_sec_tri_18, double ***Z_sec_tri_18, double ***Rzeta_sec_tri_18, double ***Zzeta_sec_tri_18, double **triangles_uv, time_t time_s, \
double ***M1, \
double *****MG, double ****M2_integral4, double ***BG, double error, \
int num_edges, int **points_sequence, double **boundary_DoFs) {
    printf("divergence_cleaning...\n");
    fflush(stdout);
    // div(B) in sections
    double ***divergence_nu_sec_tri_18 = integral_divergence_nu_dudv(G, num_sections, GLs, (int) (GL_order + 0), FuncDivergence, \
    R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, triangles_uv, num_modes, nfp_ft, B_3_tri_18N);
    timing(time_s);
    double ***divergence_sec_tri_18 = F_sec_tri_18(M1, divergence_nu_sec_tri_18, (int) num_sections, (int) triangles_uv[0][6]);
    timing(time_s);
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            delete [] divergence_nu_sec_tri_18[i][j];
        }
        delete [] divergence_nu_sec_tri_18[i];
    }
    delete [] divergence_nu_sec_tri_18;    

    // FT
    double dzeta = 2 * pi / num_sections;
    double *zeta = new double [num_sections]();
    for (int i = 0; i < num_sections; i++) {
        zeta[i] = dzeta * i;
    }

    double **divergence_tri_18N = new double *[(int) triangles_uv[0][6]];
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        divergence_tri_18N[i] = new double [18*num_modes]();
    }
    #pragma omp parallel for
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int l = 0; l < num_sections; l++) {                            
                #pragma omp atomic
                divergence_tri_18N[i][j] += divergence_sec_tri_18[l][i][j/num_modes] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
            }
        }
    }
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            delete [] divergence_sec_tri_18[i][j];
        }
        delete [] divergence_sec_tri_18[i];
    }
    delete [] divergence_sec_tri_18;    
    timing(time_s);      

    // B matrix for div(B)
    double **Bdivergence_tri_18N = new double *[(int) triangles_uv[0][6]];
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        Bdivergence_tri_18N[i] = new double [(18*num_modes)]();
    }
    #pragma omp parallel for
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        for (int k = 0; k < 18*num_modes; k++) {
            for (int t = 0; t < 18*num_modes; t++) {
                #pragma omp atomic
                Bdivergence_tri_18N[i][k] += divergence_tri_18N[i][t] * integral4(nfp_ft, num_modes, t, k, 0, 0, error) * BG[i][t/num_modes][k/num_modes];
            }   
        }
    }
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        delete [] divergence_tri_18N[i];
    }
    delete [] divergence_tri_18N;
    timing(time_s);

    // zero-valued neumann condition
    double **boundary_poisson = new double *[1];
    boundary_poisson[0] = new double [num_edges*num_modes]();

    // solve poisson equation
    double ***poisoon_1_tri_18N = solve_magnetic_field_sparse_lagrange_multiplier(Mpoisson_tri_18N_18N, Bdivergence_tri_18N, points_sequence, boundary_poisson, num_modes, (int) triangles_uv[0][6], (int) 1, (int) 1, error, time_s);
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        delete [] Bdivergence_tri_18N[i];
    }
    delete [] Bdivergence_tri_18N;
    delete [] boundary_poisson[0];
    delete [] boundary_poisson;
    timing(time_s);
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < 18*num_modes; j++) {
    //         printf("%0.9f ", poisoon_1_tri_18N[0][i][j]);
    //         if ((j+1) % (6*num_modes) == 0) {
    //             printf("    ");
    //         }
    //     }
    //     printf("\n");
    // }
    // printf("\n\n\n");
    // print_B(poisoon_1_tri_18N, G, triangles_uv, num_modes, nfp_ft, 1);
    // exit(0);
    

    // projection
    double ****projection_nu_sec_tri_18_3 = integral_projection_nu_dudv(G, num_sections, GLs, (int) (GL_order + 0), FuncProjection, \
    R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, triangles_uv, num_modes, nfp_ft, B_3_tri_18N, poisoon_1_tri_18N);
    timing(time_s);
    double ****projection_sec_tri_18_3 = F_sec_tri_18_num(M1, projection_nu_sec_tri_18_3, (int) num_sections, (int) triangles_uv[0][6], 3);
    timing(time_s);
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            for (int k = 0; k < 18; k++) {
                delete [] projection_nu_sec_tri_18_3[i][j][k];
            }
            delete [] projection_nu_sec_tri_18_3[i][j];
        }
        delete [] projection_nu_sec_tri_18_3[i];
    }
    delete [] projection_nu_sec_tri_18_3; 
    // printf("a1...\n");
    // fflush(stdout);
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
    //         for (int k = 0; k < 18; k++) {
    //             printf("%0.9f ", projection_sec_tri_18_3[0][j][k][i]);
    //             if ((k+1)%6 == 0) {
    //                 printf(" ");
    //             }
    //         }
    //         printf("\n");
    //     }
    //     printf("\n\n\n\n");
    // }   
    // printf("a2...\n");
    // fflush(stdout);
    // exit(0);

    // FT
    double ***projection_3_tri_18N = new double **[3];
    for (int i = 0; i < 3; i++) {
        projection_3_tri_18N[i] = new double *[(int) triangles_uv[0][6]];
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            projection_3_tri_18N[i][j] = new double [18*num_modes]();
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < num_sections; l++) {                            
                    #pragma omp atomic
                    projection_3_tri_18N[k][i][j] += projection_sec_tri_18_3[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
                }
            }
            
        }
    }
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            for (int k = 0; k < 18; k++) {
                delete [] projection_sec_tri_18_3[i][j][k];
            }
            delete [] projection_sec_tri_18_3[i][j];
        }
        delete [] projection_sec_tri_18_3[i];
    }
    delete [] projection_sec_tri_18_3;    
    timing(time_s);       
    // // impose dirichelet boundary conditon
    // int index1;
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         if (points_sequence[i][0+3*j] == 1) {
    //             index1 = points_sequence[i][1+3*j];
    //             for (int s = 0; s < num_modes; s++) {
    //                 projection_3_tri_18N[0][i][j * (6*num_modes) + s] = boundary_DoFs[0][index1 * (num_modes) + s];
    //                 projection_3_tri_18N[1][i][j * (6*num_modes) + s] = boundary_DoFs[1][index1 * (num_modes) + s];
    //                 projection_3_tri_18N[2][i][j * (6*num_modes) + s] = boundary_DoFs[2][index1 * (num_modes) + s];
    //             }
    //         }
    //     }
    // }  


    // // projection2
    // double ******projection_nu_sec_tri_18_2_3_3 = integral_projection2_nu_dudv(G, num_sections, GLs, (int) (GL_order + 0), FuncProjection2, \
    // R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, triangles_uv, num_modes, nfp_ft, B_3_tri_18N, poisoon_1_tri_18N);
    // timing(time_s);
    // double ******projection_sec_tri_18_2_3_3 = F_sec_tri_18_derivates_m_m(M1, projection_nu_sec_tri_18_2_3_3, (int) num_sections, (int) triangles_uv[0][6], 3, 2);
    // timing(time_s);
    // for (int i = 0; i < num_sections; i++) {
    //     for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
    //         for (int k = 0; k < 18; k++) {
    //             for (int s = 0; s < 2; s++) {
    //                 for (int t = 0; t < 3; t++) {
    //                     delete [] projection_nu_sec_tri_18_2_3_3[i][j][k][s][t];
    //                 }
    //                 delete [] projection_nu_sec_tri_18_2_3_3[i][j][k][s];
    //             }
    //             delete [] projection_nu_sec_tri_18_2_3_3[i][j][k];
    //         }
    //         delete [] projection_nu_sec_tri_18_2_3_3[i][j];
    //     }
    //     delete [] projection_nu_sec_tri_18_2_3_3[i];
    // }
    // delete [] projection_nu_sec_tri_18_2_3_3; 
    // // FT
    // double *****projection_tri_18N_2_3_3 = new double ****[(int) triangles_uv[0][6]];
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     projection_tri_18N_2_3_3[i] = new double ***[18*num_modes];
    //     for (int j = 0; j < 18*num_modes; j++) {
    //         projection_tri_18N_2_3_3[i][j] = new double **[2];
    //         for (int k = 0; k < 2; k++) {
    //             projection_tri_18N_2_3_3[i][j][k] = new double *[3];
    //             for (int l = 0; l < 3; l++) {
    //                 projection_tri_18N_2_3_3[i][j][k][l] = new double [3]();
    //             }
    //         }
    //     }
    // }
    // #pragma omp parallel for
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < 18*num_modes; j++) {
    //         for (int k = 0; k < 2; k++) {
    //             for (int s = 0; s < 3; s++) {
    //                 for (int t = 0; t < 3; t++) {  
    //                     for (int l = 0; l < num_sections; l++) {      
    //                         #pragma omp atomic
    //                         projection_tri_18N_2_3_3[i][j][k][s][t] += projection_sec_tri_18_2_3_3[l][i][j/num_modes][k][s][t] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
    //                     }
    //                 }
    //             }
    //         }
            
    //     }
    // }
    // timing(time_s);
    // for (int i = 0; i < num_sections; i++) {
    //     for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
    //         for (int k = 0; k < 18; k++) {
    //             for (int s = 0; s < 2; s++) {
    //                 for (int t = 0; t < 3; t++) {
    //                     delete [] projection_sec_tri_18_2_3_3[i][j][k][s][t];
    //                 }
    //                 delete [] projection_sec_tri_18_2_3_3[i][j][k][s];
    //             }
    //             delete [] projection_sec_tri_18_2_3_3[i][j][k];
    //         }
    //         delete [] projection_sec_tri_18_2_3_3[i][j];
    //     }
    //     delete [] projection_sec_tri_18_2_3_3[i];
    // }
    // delete [] projection_sec_tri_18_2_3_3;   
    // timing(time_s);  
    // // impose dirichelet boundary conditon
    // double ***M2_projection = M_tri_18Nn_18Nn(MG, projection_tri_18N_2_3_3, M2_integral4, (int) triangles_uv[0][6], num_modes, nfp_ft, (int) 3, error, (int) 1);
    // double **B2_projection = B_tri_18Nn(BG, projection_tri_18N_2_3_3, (int) triangles_uv[0][6], num_modes, nfp_ft, (int) 3, error, (int) 1);
    // // double ***projection_3_tri_18N = solve_magnetic_field_sparse_direct(M2_projection, B2_projection, points_sequence, boundary_DoFs, num_modes, (int) triangles_uv[0][6], (int) 3, error, time_s);
    // double ***projection_3_tri_18N = solve_magnetic_field_sparse_lagrange_multiplier(M2_projection, B2_projection, points_sequence, boundary_DoFs, num_modes, (int) triangles_uv[0][6], (int) 3, (int) 2, error, time_s);
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < (18*num_modes)*3; j++) {
    //         delete [] M2_projection[i][j];
    //     }
    //     delete [] M2_projection[i];
    //     delete [] B2_projection[i];
    // }
    // delete [] M2_projection;
    // delete [] B2_projection;

    return projection_3_tri_18N;
}


