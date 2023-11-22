#include "input_output.cpp"
#include "triangulation.cpp"
#include "refinement.cpp"
#include "meshes.cpp"
#include "integrals.cpp"
#include "integrals2.cpp"
#include "infrastructure.cpp"
#include "superstructure.cpp"
#include "solvers.cpp"
#include "solvers_SIME.cpp"
#include "header.hpp"

int main () {
    // std::ios::sync_with_stdio(false);

    // parameters
    const char *file_boundary = "input_boundary.txt";
    const char *file_coils = "coil_cfqs.txt";
    int num_sections = 120;
    int num_segments = 20;
    int nfp = 2;
    double error = 1.0E-14;
    int num_modes = 11;
    int nfp_ft = 2;
    int coil_interpolation = 29;
    int GL_order = 1; // GL_order += 1 for divergence cleaning 
    double beltrami_mu = 0.0; // 0.1 - iota = 1/3
    vector<double> initial_pressure = {20000.0, -20000.0};  // 10000.0, -10000.0

    if (num_modes % 2 == 0) {
        printf("mode number must be odd\n");
        exit(0);
    }
    time_t time_s = time(NULL);
    table();
    mn_order();
    omp_set_num_threads(NT);
    Eigen::initParallel();
    // Eigen::setNbThreads(NT);
    printf("Eigen::nbThreads = %d\n",Eigen::nbThreads());
    timing(time_s);

    
    
    // fourier coefficients
    double **fourier_coefficients = read_boundary(file_boundary);
    double **fc_zeta = fourier_coefficients_dzeta(fourier_coefficients, (int) nfp);
    double **fc_theta = fourier_coefficients_dtheta(fourier_coefficients, (int) nfp);
    double **fc_zetazeta= fourier_coefficients_dzeta(fc_zeta, (int) nfp);
    timing(time_s);
    // coils
    double ***coils_p = read_coils_p(file_coils, (int) coil_interpolation, 0.5);
    timing(time_s);
    
    // double angle = 34.0/180.0*pi;  // 90_0-30.48 91_0-28.90 20.7
    // double **triangles_RZ = ruppert_refinement((int) 80, (int) 0, (double) 0.05, 1/(2*sin(angle)), error);
    double **triangles_xy = chew_first_refinement((int) num_segments, 2.0*sin(pi/num_segments), sqrt(3), error);
    // printf("num_segments = %d, p = %d\n",num_segments, (int) triangles_xy[0][7]);
    timing(time_s);
    double ***triangles_RZ = triangles_mapping(triangles_xy, fourier_coefficients, (int) num_sections, (int) nfp);
    timing(time_s);
    triangles_reorder(triangles_RZ, triangles_xy, (int) num_sections);
    timing(time_s);
    int **points_sequence = point_retrieval(triangles_xy, num_segments, (int) triangles_xy[0][7], error);
    timing(time_s);
    double **triangles_uv = triangles_a_b_c_theta_x0_y0(triangles_xy);
    timing(time_s);
    double **boundary_DoFs = boundary_DoFs_beltrami(triangles_RZ, points_sequence, coils_p, (int) num_sections, (int) num_modes, (int) nfp_ft);
    timing(time_s);
    double **normal_3_fixed = normal_vector(fc_zeta, fc_theta, triangles_xy, points_sequence, num_sections, num_modes, nfp, nfp_ft);
    timing(time_s);

    int num_triangles = (int) triangles_uv[0][6];
    int num_vertices = points_sequence[0][10];



    double ***G = G_tri_20_18(triangles_uv);
    timing(time_s);
    double ***M1 = M_tri_18_18(G, triangles_uv);
    timing(time_s);
    double **GLs = Gauss_Legendre_quadrature_5();
    timing(time_s);

    // sparse matrix on poloidal sections
    printf("sparse matrix on poloidal section\n");

    SparseMatrix<double> M1_sparse(num_vertices * 6, num_vertices * 6);
    M1_sparse.setZero();

    double ***B_sec_tri_18 = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        B_sec_tri_18[i] = new double *[num_triangles];
        for (int j = 0; j < num_triangles; j++) {
            B_sec_tri_18[i][j] = new double [18]();
        }
    }

    typedef Eigen::Triplet<double> T;
    vector<T> coeff_M1_sparse;

    int index_j, index_k, index_J, index_K;

    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 3; j++) {
            index_j = points_sequence[i][2+3*j] * 6;
            index_J = j * 6;
            for (int k = 0; k < 3; k++) {
                index_k = points_sequence[i][2+3*k] * 6;
                index_K = k * 6;
                for (int s = 0; s < 6; s++) {
                    for (int t = 0; t < 6; t++) {
                        coeff_M1_sparse.push_back( T(index_j + s, index_k + t, M1[i][index_J + s][index_K + t]) );
                    }
                }
            }
        }
    }

    M1_sparse.setFromTriplets(coeff_M1_sparse.begin(), coeff_M1_sparse.end());
    M1_sparse.makeCompressed();
    vector<T>().swap(coeff_M1_sparse);

    solver_LU.compute(M1_sparse);
    if (solver_LU.info() != Eigen::Success) {
        printf("decomposition failed...\n");
        exit(0);
    }
    timing(time_s);
    //

    // geometry
    // R  Z  {Rzeta Zzeta  Rzetazeta Zzetazeta}
    double ****func_nu_sec_tri_18_num = integral_RZ_nu_dudv(G, num_sections, GLs, (int) GL_order, coordinate_funRZ, \
                                        (int) nfp, fourier_coefficients, fc_zeta, fc_zetazeta, (int) 2, triangles_uv);
    timing(time_s);
    // R
    #pragma omp parallel for
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            for (int k = 0; k < 18; k++) {
                B_sec_tri_18[i][j][k] = func_nu_sec_tri_18_num[i][j][k][0];
            }
        }
    }
    double ***R_sec_tri_18 = solve_surfaces_sparse(B_sec_tri_18, points_sequence, num_triangles, num_sections);
    // Z
    #pragma omp parallel for
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            for (int k = 0; k < 18; k++) {
                B_sec_tri_18[i][j][k] = func_nu_sec_tri_18_num[i][j][k][1];
            }
        }
    }
    double ***Z_sec_tri_18 = solve_surfaces_sparse(B_sec_tri_18, points_sequence, num_triangles, num_sections);    
    for (int i = 0; i < num_sections; i++) {
         for (int j = 0; j < (int) num_triangles; j++) {
            for (int k = 0; k < 18; k++) {
                delete [] func_nu_sec_tri_18_num[i][j][k];
            }
            delete [] func_nu_sec_tri_18_num[i][j];
        }
        delete [] func_nu_sec_tri_18_num[i];
    }
    delete [] func_nu_sec_tri_18_num;

    // for (int i = 0; i < num_triangles; i++) {
    //    printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n", \
    //    R_sec_tri_18[40][i][0],Z_sec_tri_18[40][i][0], R_sec_tri_18[40][i][6],Z_sec_tri_18[40][i][6], R_sec_tri_18[40][i][12],Z_sec_tri_18[40][i][12]);
    // }
    // exit(0);

    // Rzeta Zzeta Rzetazeta Zzetazeta
    double ***Rzeta_sec_tri_18 = new double **[num_sections];
    double ***Zzeta_sec_tri_18 = new double **[num_sections];
    double ***Rzetazeta_sec_tri_18 = new double **[num_sections];
    double ***Zzetazeta_sec_tri_18 = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        Rzeta_sec_tri_18[i] = new double *[num_triangles];
        Zzeta_sec_tri_18[i] = new double *[num_triangles];
        Rzetazeta_sec_tri_18[i] = new double *[num_triangles];
        Zzetazeta_sec_tri_18[i] = new double *[num_triangles];
        for (int j = 0; j < num_triangles; j++) {
            Rzeta_sec_tri_18[i][j] = new double [18]();
            Zzeta_sec_tri_18[i][j] = new double [18]();
            Rzetazeta_sec_tri_18[i][j] = new double [18]();
            Zzetazeta_sec_tri_18[i][j] = new double [18]();
        }
    }

    double dzeta = 2 * pi / num_sections;
    double *zeta = new double [num_sections]();
    for (int i = 0; i < num_sections; i++) {
        zeta[i] = dzeta * i;
    }

    // dzeta d^2zeta;
    double **R_tri_18N = new double *[num_triangles];
    double **Z_tri_18N = new double *[num_triangles];
    for (int i = 0; i < num_triangles; i++) {
        R_tri_18N[i] = new double [18*num_modes]();
        Z_tri_18N[i] = new double [18*num_modes]();
    }
    #pragma omp parallel for
    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < num_sections; k++) {
                #pragma omp atomic
                R_tri_18N[i][j] += R_sec_tri_18[k][i][j/num_modes] * fourier_basis(j%num_modes, (int) nfp, zeta[k]) * dzeta;
                #pragma omp atomic
                Z_tri_18N[i][j] += Z_sec_tri_18[k][i][j/num_modes] * fourier_basis(j%num_modes, (int) nfp, zeta[k]) * dzeta;
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            for (int k = 0; k < 18*num_modes; k++) {
                #pragma omp atomic
                Rzeta_sec_tri_18[i][j][k/num_modes] += R_tri_18N[j][k] * fourier_basis_test2(k%num_modes, nfp, zeta[i]);
                #pragma omp atomic
                Zzeta_sec_tri_18[i][j][k/num_modes] += Z_tri_18N[j][k] * fourier_basis_test2(k%num_modes, nfp, zeta[i]);
                #pragma omp atomic
                Rzetazeta_sec_tri_18[i][j][k/num_modes] += R_tri_18N[j][k] * fourier_basis_test3(k%num_modes, nfp, zeta[i]);
                #pragma omp atomic
                Zzetazeta_sec_tri_18[i][j][k/num_modes] += Z_tri_18N[j][k] * fourier_basis_test3(k%num_modes, nfp, zeta[i]);
            }
        }
    }
    timing(time_s);
    // for (int i = 0; i < num_triangles; i++) {
    //     for (int j = 0; j < 18; j++) {
    //         printf("%0.9f ",Zzeta_sec_tri_18[0][i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n\n\n\n");
    // for (int i = 0; i < num_triangles; i++) {
    //     for (int j = 0; j < 18; j++) {
    //         printf("%0.9f ",Zzeta_sec_tri_18[40][i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n\n\n\n");
    // for (int i = 0; i < num_triangles; i++) {
    //     for (int j = 0; j < 18; j++) {
    //         printf("%0.9f ",Zzetazeta_sec_tri_18[0][i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n\n\n\n");
    // for (int i = 0; i < num_triangles; i++) {
    //     for (int j = 0; j < 18; j++) {
    //         printf("%0.9f ",Zzetazeta_sec_tri_18[40][i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n\n\n\n");
    // exit(0);
    //
    
    // divergence constriant
    double ****divcoeff_sec_tri_3_var = divergence_free_constraints(triangles_uv, num_sections, num_modes, nfp_ft, R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18);
    timing(time_s);
    //

    // pressure initialization
    double ***pressure_nu_sec_tri_18 = integral_pressure_nu_dudv(G, num_sections, GLs, (int) GL_order, pressure_initialization, triangles_uv, initial_pressure);
    timing(time_s);

    double ***pressure_sec_tri_18 = solve_surfaces_sparse(pressure_nu_sec_tri_18, points_sequence, num_triangles, num_sections);   
    for (int i = 0; i < num_sections; i++) {
         for (int j = 0; j < (int) num_triangles; j++) {
            delete [] pressure_nu_sec_tri_18[i][j];          
        }
        delete [] pressure_nu_sec_tri_18[i];
    }
    delete [] pressure_nu_sec_tri_18;
    //

    // jacobian
    double ***J_nu_sec_tri_18 = integral_jacobian_nu_dudv(G, num_sections, GLs, (int) GL_order, jacobian, R_sec_tri_18, Z_sec_tri_18, triangles_uv);
    timing(time_s);

    double ***J_sec_tri_18 = solve_surfaces_sparse(J_nu_sec_tri_18, points_sequence, num_triangles, num_sections);
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            delete [] J_nu_sec_tri_18[i][j];
        }
        delete [] J_nu_sec_tri_18[i];
    }
    delete [] J_nu_sec_tri_18;
    //

    // beltrami field
    int num_scalar = 3;
    int num_operators = 11;
    double ***temp_sec_tri_18;
    double ******beltrami_nu_sec_tri_18_11_sca_sca = integral_beltrami_nu_dudv(G, num_sections, GLs, (int) GL_order, FuncBeltrami, \
    R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, Rzetazeta_sec_tri_18, Zzetazeta_sec_tri_18, triangles_uv, (int) num_scalar, beltrami_mu);
    timing(time_s);

    double ******beltrami_sec_tri_18_11_sca_sca = new double *****[num_sections];
    for (int i = 0; i < num_sections; i++) {
        beltrami_sec_tri_18_11_sca_sca[i] = new double ****[num_triangles];
        for (int j = 0; j < num_triangles; j++) {
            beltrami_sec_tri_18_11_sca_sca[i][j] = new double ***[18];
            for (int k = 0; k < 18; k++) {
                beltrami_sec_tri_18_11_sca_sca[i][j][k] = new double **[num_operators];
                for (int s = 0; s < num_operators; s++) {
                    beltrami_sec_tri_18_11_sca_sca[i][j][k][s] = new double *[num_scalar];
                    for (int t = 0; t < num_scalar; t++) {
                        beltrami_sec_tri_18_11_sca_sca[i][j][k][s][t] = new double [num_scalar]();
                    }
                }
            }
        }
    }

    for (int i = 0; i < num_operators; i++) {
        for (int j = 0; j < num_scalar; j++) {
            for (int k = 0; k < num_scalar; k++) {
                
                #pragma omp parallel for
                for (int l = 0; l < num_sections; l++) {
                    for (int s = 0; s < num_triangles; s++) {
                        for (int t = 0; t < 18; t++) {
                            B_sec_tri_18[l][s][t] = beltrami_nu_sec_tri_18_11_sca_sca[l][s][t][i][j][k];
                        }
                    }
                }

                temp_sec_tri_18 = solve_surfaces_sparse(B_sec_tri_18, points_sequence, num_triangles, num_sections);

                #pragma omp parallel for
                for (int l = 0; l < num_sections; l++) {
                    for (int s = 0; s < num_triangles; s++) {
                        for (int t = 0; t < 18; t++) {
                            beltrami_sec_tri_18_11_sca_sca[l][s][t][i][j][k] = temp_sec_tri_18[l][s][t];
                        }
                    }
                }

                for (int l = 0; l < num_sections; l++) {
                    for (int s = 0; s < num_triangles; s++) {
                        delete [] temp_sec_tri_18[l][s];
                    }
                    delete [] temp_sec_tri_18[l];
                }
                delete [] temp_sec_tri_18;


            }
        }
    }
    timing(time_s);

    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            for (int k = 0; k < 18; k++) {
                for (int s = 0; s < 11; s++) {
                    for (int t = 0; t < num_scalar; t++) {
                        delete [] beltrami_nu_sec_tri_18_11_sca_sca[i][j][k][s][t];
                    }
                    delete [] beltrami_nu_sec_tri_18_11_sca_sca[i][j][k][s];
                }
                delete [] beltrami_nu_sec_tri_18_11_sca_sca[i][j][k];
            }
            delete [] beltrami_nu_sec_tri_18_11_sca_sca[i][j];
        }
        delete [] beltrami_nu_sec_tri_18_11_sca_sca[i];
    }
    delete [] beltrami_nu_sec_tri_18_11_sca_sca;
    //

    // poisson equation
    int num_poisson = 10;
    double ****poisson_nu_sec_tri_18_10 = integral_poisson_nu_dudv(G, num_sections, GLs, (int) GL_order, FuncPoisson, \
    R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, Rzetazeta_sec_tri_18, Zzetazeta_sec_tri_18, triangles_uv, (int) 10);
    timing(time_s);

    double ****poisson_sec_tri_18_10 = solve_surfaces_sparse_num(poisson_nu_sec_tri_18_10, temp_sec_tri_18, B_sec_tri_18, points_sequence, num_sections, num_triangles, (int) num_poisson);
    timing(time_s);

    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            for (int k = 0; k < 18; k++) {
                delete [] poisson_nu_sec_tri_18_10[i][j][k];
            }
            delete [] poisson_nu_sec_tri_18_10[i][j];
        }
        delete [] poisson_nu_sec_tri_18_10[i];
    }
    delete [] poisson_nu_sec_tri_18_10;
    //

    // double ****divcoeff_nu_sec_tri_18_9 = integral_DivCoeff_nu_dudv(G, num_sections, GLs, (int) GL_order, FuncDivCoeff, \
    // R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, triangles_uv);
    // timing(time_s);
    // double ****divcoeff_sec_tri_18_9 = F_sec_tri_18_num(M1, divcoeff_nu_sec_tri_18_9, (int) num_sections, (int) triangles_uv[0][6], (int) 9);
    // timing(time_s);
    // for (int i = 0; i < num_sections; i++) {
    //     for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
    //         for (int k = 0; k < 18; k++) {
    //             delete [] divcoeff_nu_sec_tri_18_9[i][j][k];
    //         }
    //         delete [] divcoeff_nu_sec_tri_18_9[i][j];
    //     }
    //     delete [] divcoeff_nu_sec_tri_18_9[i];
    // }
    // delete [] divcoeff_nu_sec_tri_18_9; 


    // // test for FT
    // int num_m = num_sections / 4;
    // double dzzeta = 2 * pi / num_sections, zzeta = 0;
    // double **R_tri_18N = new double *[(int) triangles_uv[0][6]];
    // double **Z_tri_18N = new double *[(int) triangles_uv[0][6]];
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     R_tri_18N[i] = new double [18*num_m];
    //     Z_tri_18N[i] = new double [18*num_m];
    //     for (int j = 0; j < 18*num_m; j++) {
    //         R_tri_18N[i][j] = 0;
    //         Z_tri_18N[i][j] = 0;
    //         for (int k = 0; k < num_sections; k++) {
    //             zzeta = dzzeta * k;
    //             R_tri_18N[i][j] += R_sec_tri_18[k][i][j/num_m] * fourier_basis(j%num_m, (int) 1, zzeta) *dzzeta;
    //             Z_tri_18N[i][j] += Z_sec_tri_18[k][i][j/num_m] * fourier_basis(j%num_m, (int) 1, zzeta) *dzzeta;
    //         }
    //     }
    // }
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < 18*num_m; j++) {
    //         printf("%0.9f  ", R_tri_18N[i][j]);
    //         if (j%num_m == num_m-1) {
    //             printf("\n");
    //         }
    //     }
    //     printf("\n");
    // }
    // printf("\n\n\n\n");
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < 18*num_m; j++) {
    //         printf("%0.9f  ", Z_tri_18N[i][j]);
    //         if (j%num_m == num_m-1) {
    //             printf("\n");
    //         }
    //     }
    //     printf("\n");
    // }
    // printf("\n\n\n\n");
    // //

    /*
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            delete [] R_sec_tri_18[i][j];      
            delete [] Z_sec_tri_18[i][j]; 
            delete [] Rzeta_sec_tri_18[i][j]; 
            delete [] Zzeta_sec_tri_18[i][j]; 
            delete [] Rzetazeta_sec_tri_18[i][j]; 
            delete [] Zzetazeta_sec_tri_18[i][j];       
        }
        delete [] R_sec_tri_18[i];
        delete [] Z_sec_tri_18[i];
        delete [] Rzeta_sec_tri_18[i];
        delete [] Zzeta_sec_tri_18[i];
        delete [] Rzetazeta_sec_tri_18[i];
        delete [] Zzetazeta_sec_tri_18[i];
    }
    delete [] R_sec_tri_18;
    delete [] Z_sec_tri_18;
    delete [] Rzeta_sec_tri_18;
    delete [] Zzeta_sec_tri_18;
    delete [] Rzetazeta_sec_tri_18;
    delete [] Zzetazeta_sec_tri_18;
    */

    // fourier transform
    printf("FT...\n");

    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < 18*num_modes; j++) {
    //         if (j % (6*num_modes) == 0) {
    //             printf("\n");
    //         }
    //         printf("%0.9f ", func_tri_18N_num[i][j][0]);
    //     }
    //     printf("\n");
    // }
    // printf("1...\n");
    // fflush(stdout);
    // exit(0);

    // // test for fourier transfrom
    // double uu[4], vv[4], uv[4][20], rr[4], zz[4], r_modes[18], z_modes[18], zzeta = pi / 2;
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < 18; j++) {
    //         r_modes[j] = 0;
    //         z_modes[j] = 0;
    //     }
    //     for (int j = 0; j < 18*num_modes; j++) {
    //         r_modes[j/num_modes] += func_tri_18N_num[i][j][0] * fourier_basis_test(j%num_modes, (int) nfp_ft, zzeta);
    //         z_modes[j/num_modes] += func_tri_18N_num[i][j][1] * fourier_basis_test(j%num_modes, (int) nfp_ft, zzeta);
    //     }
    //     uu[0] = (triangles_uv[i][0]-triangles_uv[i][1]) / 3.0; // barycenter
    //     vv[0] = (triangles_uv[i][2]) / 3.0;
    //     uu[1] = -triangles_uv[i][1];
    //     vv[1] = 0;
    //     uu[2] = triangles_uv[i][0]; 
    //     vv[2] = 0;
    //     uu[3] = 0; 
    //     vv[3] = triangles_uv[i][2];
    //     for (int l = 0; l < 4; l++) {
    //         rr[l] = 0;
    //         zz[l] = 0;
    //         for (int j = 0; j < 20; j++) {
    //             uv[l][j] = power_nonnegative(uu[l], m_order[j]) *   power_nonnegative(vv[l], n_order[j]);
    //             for (int k = 0; k < 18; k++) {
    //                 rr[l] += G[i][j][k] * uv[l][j] * r_modes[k];
    //                 zz[l] += G[i][j][k] * uv[l][j] * z_modes[k];
    //             }
    //         }
    //         printf("%0.9f  %0.9f\n", rr[l], zz[l]);
    //     }  
    // }
    // exit(0);

    double **pressure_tri_18N = new double *[num_triangles];
    for (int i = 0; i < num_triangles; i++) {
        pressure_tri_18N[i] = new double [18*num_modes]();
        for (int j = 0; j < 18*num_modes; j++) {
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < num_sections; k++) {
                #pragma omp atomic
                pressure_tri_18N[i][j] += pressure_sec_tri_18[k][i][j/num_modes] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[k]) * dzeta;
            }
        }
    }
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            delete [] pressure_sec_tri_18[i][j];
        }
        delete [] pressure_sec_tri_18[i];
    }
    delete [] pressure_sec_tri_18;
    timing(time_s); 
    // for (int i = 0; i < num_triangles; i++) {
    //     for (int j = 0; j < 18; j++) {
    //         printf("%0.9f ", pressure_tri_18N[i][j*num_modes]);
    //         if ((j+1)%6 == 0) {
    //             printf("     ");
    //         }
    //     }
    //     printf("\n");
    // }
    // exit(0);
    // print_p(pressure_tri_18N, G, triangles_uv, num_modes, nfp_ft);
    // exit(0);

    double **J_tri_18N = new double *[num_triangles];
    for (int i = 0; i < num_triangles; i++) {
        J_tri_18N[i] = new double [18*num_modes]();
        for (int j = 0; j < 18*num_modes; j++) {
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < num_sections; k++) {
                #pragma omp atomic
                J_tri_18N[i][j] += J_sec_tri_18[k][i][j/num_modes] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[k]) * dzeta;
            }
        }
    }
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            delete [] J_sec_tri_18[i][j];
        }
        delete [] J_sec_tri_18[i];
    }
    delete [] J_sec_tri_18;
    timing(time_s);  


    double *****beltrami_tri_18N_11_sca_sca = new double ****[num_triangles];
    for (int i = 0; i < num_triangles; i++) {
        beltrami_tri_18N_11_sca_sca[i] = new double ***[18*num_modes];
        for (int j = 0; j < 18*num_modes; j++) {
            beltrami_tri_18N_11_sca_sca[i][j] = new double **[11];
            for (int k = 0; k < 11; k++) {
                beltrami_tri_18N_11_sca_sca[i][j][k] = new double *[num_scalar];
                for (int s = 0; s < num_scalar; s++) {
                    beltrami_tri_18N_11_sca_sca[i][j][k][s] = new double [num_scalar]();
                }
            }
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < 11; k++) {
                for (int s = 0; s < num_scalar; s++) {
                    for (int t = 0; t < num_scalar; t++) {
                        for (int l = 0; l < num_sections; l++) {                            
                            #pragma omp atomic
                            beltrami_tri_18N_11_sca_sca[i][j][k][s][t] += beltrami_sec_tri_18_11_sca_sca[l][i][j/num_modes][k][s][t] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            for (int k = 0; k < 18; k++) {
                for (int s = 0; s < 11; s++) {
                    for (int t = 0; t < num_scalar; t++) {
                        delete [] beltrami_sec_tri_18_11_sca_sca[i][j][k][s][t];
                    }
                    delete [] beltrami_sec_tri_18_11_sca_sca[i][j][k][s];
                }
                delete [] beltrami_sec_tri_18_11_sca_sca[i][j][k];
            }
            delete [] beltrami_sec_tri_18_11_sca_sca[i][j];
        }
        delete [] beltrami_sec_tri_18_11_sca_sca[i];
    }
    delete [] beltrami_sec_tri_18_11_sca_sca;    
    timing(time_s);     


    double ***poisson_tri_18N_10 = new double **[num_triangles];
    for (int i = 0; i < num_triangles; i++) {
        poisson_tri_18N_10[i] = new double *[18*num_modes];
        for (int j = 0; j < 18*num_modes; j++) {
            poisson_tri_18N_10[i][j] = new double [10]();
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < 10; k++) {
                for (int l = 0; l < num_sections; l++) {                            
                    #pragma omp atomic
                    poisson_tri_18N_10[i][j][k] += poisson_sec_tri_18_10[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
                }
            }
        }
    }
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            for (int k = 0; k < 18; k++) {
                delete [] poisson_sec_tri_18_10[i][j][k];
            }
            delete [] poisson_sec_tri_18_10[i][j];
        }
        delete [] poisson_sec_tri_18_10[i];
    }
    delete [] poisson_sec_tri_18_10;    
    timing(time_s);    

    // double ****divcoeff_tri_3_var_modes = new double ***[num_triangles];
    // for (int i = 0; i < num_triangles; i++) {
    //     divcoeff_tri_3_var_modes[i] = new double **[3];
    //     for (int j = 0; j < 3; j++) {
    //         divcoeff_tri_3_var_modes[i][j] = new double *[9];
    //         for (int k = 0; k < 9; k++) {
    //             divcoeff_tri_3_var_modes[i][j][k] = new double [num_modes]();
    //         }
    //     }
    // }
    // #pragma omp parallel for
    // for (int i = 0; i < num_triangles; i++) {
    //     for (int j = 0; j < 18*num_modes; j++) {
    //         for (int k = 0; k < 9; k++) {
    //             for (int l = 0; l < num_sections; l++) {                            
    //                 #pragma omp atomic
    //                 divcoeff_tri_3_var_modes[i][j/(6*num_modes)][k][j%num_modes] += divcoeff_sec_tri_18_9[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
    //             }
    //         }
    //     }
    // }
    // for (int i = 0; i < num_sections; i++) {
    //     for (int j = 0; j < num_triangles; j++) {
    //         for (int k = 0; k < 18; k++) {
    //             delete [] divcoeff_sec_tri_18_9[i][j][k];
    //         }
    //         delete [] divcoeff_sec_tri_18_9[i][j];
    //     }
    //     delete [] divcoeff_sec_tri_18_9[i];
    // }
    // delete [] divcoeff_sec_tri_18_9;    
    // timing(time_s);    

       
    double *****MG = MG_tri_6_18_18_18(G, triangles_uv);
    timing(time_s);
    double ***BG = BG_tri_18_18(G, triangles_uv);
    timing(time_s);
    double ****M2_integral4 = M_integral4(nfp_ft, num_modes, error);
    timing(time_s);
    double ***M2_poisson = Mpoisson_tri_18N_18N(MG, poisson_tri_18N_10, M2_integral4, num_triangles, num_modes);
    timing(time_s);
    double ***M2 = M_tri_18Nn_18Nn(MG, beltrami_tri_18N_11_sca_sca, M2_integral4, num_triangles, num_modes, nfp_ft, num_scalar, error, num_operators-1);
    timing(time_s);
    double **B2 = B_tri_18Nn(BG, beltrami_tri_18N_11_sca_sca, num_triangles, num_modes, nfp_ft, num_scalar, error, num_operators-1);
    timing(time_s);
    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < 11; k++) {
                for (int l = 0; l < num_scalar; l++) {
                    delete [] beltrami_tri_18N_11_sca_sca[i][j][k][l];
                }
                delete [] beltrami_tri_18N_11_sca_sca[i][j][k];
            }
            delete [] beltrami_tri_18N_11_sca_sca[i][j];
        }
        delete [] beltrami_tri_18N_11_sca_sca[i];
    }
    delete [] beltrami_tri_18N_11_sca_sca;   
    // exit(0);
    // solve_poisson_sparse(M2, B2, points_sequence, num_modes, num_triangles, time_s);
    // double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_direct(M2, B2, points_sequence, boundary_DoFs, num_modes, num_triangles, (int) num_scalar, error, time_s);
    // double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_lagrange_multiplier(M2, B2, points_sequence, boundary_DoFs, num_modes, num_triangles, (int) num_scalar, (int) 0, error, time_s);
    // double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_lagrange_multiplier_div(divcoeff_sec_tri_3_var, M2, B2, points_sequence, boundary_DoFs, num_modes, num_triangles, (int) nfp_ft, (int) num_scalar, (int) 0, error, time_s);
    // double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_lagrange_multiplier_div(divcoeff_sec_tri_3_var, M2, B2, points_sequence, normal_3_fixed, num_modes, num_triangles, (int) nfp_ft, (int) num_scalar, (int) 1, error, time_s);

    // double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_ALM_div(divcoeff_sec_tri_3_var, M2, B2, points_sequence, boundary_DoFs, normal_3_fixed, num_modes, num_triangles, (int) nfp_ft, (int) num_scalar, (int) 0, error, time_s);
    // double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_ALM_div(divcoeff_sec_tri_3_var, M2, B2, points_sequence, boundary_DoFs, normal_3_fixed, num_modes, num_triangles, (int) nfp_ft, (int) num_scalar, (int) 1, error, time_s);

    double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_LM(M2, B2, points_sequence, boundary_DoFs, num_modes, num_triangles, (int) num_scalar, error, time_s);

    for (int i = 0; i < num_triangles; i++) {
        delete [] B2[i];
    }
    delete [] B2;
    timing(time_s);
    // print_RZ(R_tri_18N, Z_tri_18N, G, triangles_uv, num_modes, nfp_ft);
    // print_p(pressure_tri_18N, G, triangles_uv, num_modes, nfp_ft);
    // print_B(beltrami_3_tri_18N, G, triangles_uv, num_modes, nfp_ft, num_scalar);
    // exit(0);

    double ***mf1 = beltrami_3_tri_18N;
    double ***mf2;

























    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // SIME
    printf("solving for static ideal mhd equilibrium...\n");
    double norm1_B[3], norm2_B[3];
    double norm1_gradP[3], norm2_gradP[3];
    double norm1_pres, norm2_pres2;
    double blending = 1.0;

    double **pres1 = pressure_tri_18N, ***pres2;

    double ****Jperpendicular_nu_sec_tri_18_3;
    double ****Jperpendicular_sec_tri_18_3;
    double ***Jperpendicular_3_tri_18N;

    double ****GradSigma_nu_sec_tri_18_3;
    double ****GradSigma_sec_tri_18_3;
    double ***GradSigma_3_tri_18N;

    double ***DivGradSigma_nu_sec_tri_18;
    double ***DivGradSigma_sec_tri_18;
    double **DivGradSigma_tri_18N; 

    double **B2_DivGradSigma;

    double ***sigma_1_tri_18N; 

    double ****Jparallel_nu_sec_tri_18_3;
    double ****Jparallel_sec_tri_18_3;
    double ***Jparallel_3_tri_18N;

    double ***J_3_tri_18N;
    double ***J1 = new double **[3];
    for (int i = 0; i < 3; i++) {
        J1[i] = new double *[num_triangles];
        for (int j = 0; j < num_triangles; j++) {
            J1[i][j] = new double [18*num_modes]();
        }
    }

    double ****curlJ_nu_sec_tri_18_3;
    double ****curlJ_sec_tri_18_3;
    double ***curlJ_3_tri_18N;
    double **B2_curlJ;

    // grad(pressure)_parallel
    double ****gradP_nu_sec_tri_18_3 = integral_gradP_nu_dudv(G, num_sections, num_modes, nfp_ft, GLs, (int) GL_order, FuncGradP, \
                                                              R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, pres1, triangles_uv);
    timing(time_s);
    double ****gradP_sec_tri_18_3 = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        gradP_sec_tri_18_3[i] = new double **[num_triangles];
        for (int j = 0; j < num_triangles; j++) {
            gradP_sec_tri_18_3[i][j] = new double *[18];
            for (int k = 0; k < 18; k++) {
                gradP_sec_tri_18_3[i][j][k] = new double [3]();
            }
        }
    }

    gradP_sec_tri_18_3 = solve_surfaces_sparse_num(gradP_nu_sec_tri_18_3, temp_sec_tri_18, B_sec_tri_18, points_sequence, num_sections, num_triangles, (int) 3);
    // timing(time_s);

    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            for (int k = 0; k < 18; k++) {
                delete [] gradP_nu_sec_tri_18_3[i][j][k];
            }
            delete [] gradP_nu_sec_tri_18_3[i][j];
        }
        delete [] gradP_nu_sec_tri_18_3[i];
    }
    delete [] gradP_nu_sec_tri_18_3;  
    // FT
    double ***gradP_3_tri_18N = new double **[3];
    for (int i = 0; i < 3; i++) {
        gradP_3_tri_18N[i] = new double *[num_triangles];
        for (int j = 0; j < num_triangles; j++) {
            gradP_3_tri_18N[i][j] = new double [18*num_modes]();
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < num_sections; l++) {                            
                    #pragma omp atomic
                    gradP_3_tri_18N[k][i][j] += gradP_sec_tri_18_3[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
                }
            }
        }
    }
    timing(time_s);   
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            for (int k = 0; k < 18; k++) {
                delete [] gradP_sec_tri_18_3[i][j][k];
            }
            delete [] gradP_sec_tri_18_3[i][j];
        }
        delete [] gradP_sec_tri_18_3[i];
    }
    delete [] gradP_sec_tri_18_3;
    // print_B(gradP_3_tri_18N, G, triangles_uv, num_modes, nfp_ft, num_scalar);
    //

    double ***gradP1 = gradP_3_tri_18N, ***gradP2;

    double ***divgradPparallel_nu_sec_tri_18;
    double ***divgradPparallel_sec_tri_18;
    double **divgradPparallel_tri_18N;
    double **B2_divgradPperp;

    double **boundary_zeros = new double *[3];
    for (int i = 0; i < 3; i++) {
        boundary_zeros[i] = new double [num_segments*num_modes]();
    }

    double **boundary_pressure = new double *[1];
    boundary_pressure[0] = new double [num_segments*num_modes]();
    int index_pres1;
    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 3; j++) {
            if (points_sequence[i][0+3*j] == 1) {
                index_pres1 = points_sequence[i][1+3*j]*num_modes;
                for (int k = 0; k < num_modes; k++) {
                    boundary_pressure[0][index_pres1+k] = pressure_tri_18N[i][j*(6*num_modes)+k];
                }
            }
        }
    }

    // for (int i = 0; i < num_triangles; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         for (int k = 0; k < num_modes; k++) {
    //             printf("%0.9f ",pressure_tri_18N[i][j*(6*num_modes)+k]);
    //         }
    //         printf("        ");
    //     }
    //     printf("\n");
    // }
    // printf("\n\n");

    // for (int i = 0; i < num_segments*num_modes; i++) {
    //     printf("%0.9f ",boundary_pressure[0][i]);
    //     if ((i+1) % num_modes == 0) {
    //         printf("\n");
    //     }
    // }
    // printf("\n");
    // exit(0);

    for (int iter = 0; iter < 1000; iter++) {
        // printf("SIME iteration: %d\n", iter + 1);

        // ---------------------------------------- iteration of grad(P)
        gradP_nu_sec_tri_18_3 = integral_gradPparallel_nu_dudv(G, num_sections, num_modes, nfp_ft, GLs, (int) GL_order, FuncGradPparallel, \
        R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, mf1, gradP1, triangles_uv);
        timing(time_s);
        gradP_sec_tri_18_3 = new double ***[num_sections];
        for (int i = 0; i < num_sections; i++) {
            gradP_sec_tri_18_3[i] = new double **[num_triangles];
            for (int j = 0; j < num_triangles; j++) {
                gradP_sec_tri_18_3[i][j] = new double *[18];
                for (int k = 0; k < 18; k++) {
                    gradP_sec_tri_18_3[i][j][k] = new double [3]();
                }
            }
        }

        gradP_sec_tri_18_3 = solve_surfaces_sparse_num(gradP_nu_sec_tri_18_3, temp_sec_tri_18, B_sec_tri_18, points_sequence, num_sections, num_triangles, (int) 3);
        // timing(time_s);

        for (int i = 0; i < num_sections; i++) {
            for (int j = 0; j < num_triangles; j++) {
                for (int k = 0; k < 18; k++) {
                    delete [] gradP_nu_sec_tri_18_3[i][j][k];
                }
                delete [] gradP_nu_sec_tri_18_3[i][j];
            }
            delete [] gradP_nu_sec_tri_18_3[i];
        }
        delete [] gradP_nu_sec_tri_18_3;  
        // FT
        gradP2 = new double **[3];
        for (int i = 0; i < 3; i++) {
            gradP2[i] = new double *[num_triangles];
            for (int j = 0; j < num_triangles; j++) {
                gradP2[i][j] = new double [18*num_modes]();
            }
        }
        #pragma omp parallel for
        for (int i = 0; i < num_triangles; i++) {
            for (int j = 0; j < 18*num_modes; j++) {
                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < num_sections; l++) {                            
                        #pragma omp atomic
                        gradP2[k][i][j] += gradP_sec_tri_18_3[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
                    }
                }
            }
        }
        // timing(time_s);
        for (int i = 0; i < num_sections; i++) {
            for (int j = 0; j < num_triangles; j++) {
                for (int k = 0; k < 18; k++) {
                    delete [] gradP_sec_tri_18_3[i][j][k];
                }
                delete [] gradP_sec_tri_18_3[i][j];
            }
            delete [] gradP_sec_tri_18_3[i];
        }
        delete [] gradP_sec_tri_18_3;
        // gradP.norm()
        for (int i = 0; i < 3; i++) {
            norm1_gradP[i] = 0;
            norm2_gradP[i] = 0;
            for (int j = 0; j < num_triangles; j++) {
                for (int k = 0; k < 18*num_modes; k++) {
                    norm1_gradP[i] += gradP1[i][j][k] * gradP1[i][j][k];
                    norm2_gradP[i] += gradP2[i][j][k] * gradP2[i][j][k];
                }
            }
            norm1_gradP[i] /= num_triangles * (18*num_modes);
            norm2_gradP[i] /= num_triangles * (18*num_modes);
        }
        printf("iter: %3d, gradP, estimated error = %0.6e,  %0.6e,  %0.6e\n", iter + 1, abs(norm2_gradP[0]-norm1_gradP[0]) / norm1_gradP[0], abs(norm2_gradP[1]-norm1_gradP[1]) / norm1_gradP[1], abs(norm2_gradP[2]-norm1_gradP[2]) / norm1_gradP[2]);
        fflush(stdout);
        // exchange
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < num_triangles; j++) {
                for (int k = 0; k < 18*num_modes; k++) {
                    gradP1[i][j][k] = blending * gradP2[i][j][k] + (1-blending) * gradP1[i][j][k];
                }
            }
        }

        // print_B(gradP1, G, triangles_uv, num_modes, nfp_ft, num_scalar);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < num_triangles; j++) {
                delete [] gradP2[i][j];
            }
            delete [] gradP2[i];
        }
        delete [] gradP2;

        // // ---------------------------------------- iteration of J
        // // perpendicular current density
        // Jperpendicular_nu_sec_tri_18_3 = integral_Jperpendicular_nu_dudv(G, num_sections, num_modes, nfp_ft, GLs, (int) GL_order, FuncJperpendicular, \
        // R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, mf1, gradP1, triangles_uv);
        // timing(time_s);

        // Jperpendicular_sec_tri_18_3 = solve_surfaces_sparse_num(Jperpendicular_nu_sec_tri_18_3, temp_sec_tri_18, B_sec_tri_18, points_sequence, num_sections, num_triangles, (int) 3);
        // // timing(time_s);

        // for (int i = 0; i < num_sections; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18; k++) {
        //             delete [] Jperpendicular_nu_sec_tri_18_3[i][j][k];
        //         }
        //         delete [] Jperpendicular_nu_sec_tri_18_3[i][j];
        //     }
        //     delete [] Jperpendicular_nu_sec_tri_18_3[i];
        // }
        // delete [] Jperpendicular_nu_sec_tri_18_3;  
        // // FT
        // Jperpendicular_3_tri_18N = new double **[3];
        // for (int i = 0; i < 3; i++) {
        //     Jperpendicular_3_tri_18N[i] = new double *[num_triangles];
        //     for (int j = 0; j < num_triangles; j++) {
        //         Jperpendicular_3_tri_18N[i][j] = new double [18*num_modes]();
        //     }
        // }
        // #pragma omp parallel for
        // for (int i = 0; i < num_triangles; i++) {
        //     for (int j = 0; j < 18*num_modes; j++) {
        //         for (int k = 0; k < 3; k++) {
        //             for (int l = 0; l < num_sections; l++) {                            
        //                 #pragma omp atomic
        //                 Jperpendicular_3_tri_18N[k][i][j] += Jperpendicular_sec_tri_18_3[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
        //             }
        //         }
        //     }
        // }
        // // timing(time_s);  
        // for (int i = 0; i < num_sections; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18; k++) {
        //             delete [] Jperpendicular_sec_tri_18_3[i][j][k];
        //         }
        //         delete [] Jperpendicular_sec_tri_18_3[i][j];
        //     }
        //     delete [] Jperpendicular_sec_tri_18_3[i];
        // }
        // delete [] Jperpendicular_sec_tri_18_3;  

        // // parallel current density
        // // GradSigma
        // GradSigma_nu_sec_tri_18_3 = integral_gradSigma_nu_dudv(G, num_sections, num_modes, nfp_ft, GLs, (int) GL_order, FuncGradSigma, \
        // R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, Jperpendicular_3_tri_18N, mf1, triangles_uv);
        // timing(time_s);

        // GradSigma_sec_tri_18_3 = solve_surfaces_sparse_num(GradSigma_nu_sec_tri_18_3, temp_sec_tri_18, B_sec_tri_18, points_sequence, num_sections, num_triangles, (int) 3);
        // // timing(time_s);

        // for (int i = 0; i < num_sections; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18; k++) {
        //             delete [] GradSigma_nu_sec_tri_18_3[i][j][k];
        //         }
        //         delete [] GradSigma_nu_sec_tri_18_3[i][j];
        //     }
        //     delete [] GradSigma_nu_sec_tri_18_3[i];
        // }
        // delete [] GradSigma_nu_sec_tri_18_3;  
        // // FT
        // GradSigma_3_tri_18N = new double **[3];
        // for (int i = 0; i < 3; i++) {
        //     GradSigma_3_tri_18N[i] = new double *[num_triangles];
        //     for (int j = 0; j < num_triangles; j++) {
        //         GradSigma_3_tri_18N[i][j] = new double [18*num_modes]();
        //     }
        // }
        // #pragma omp parallel for
        // for (int i = 0; i < num_triangles; i++) {
        //     for (int j = 0; j < 18*num_modes; j++) {
        //         for (int k = 0; k < 3; k++) {
        //             for (int l = 0; l < num_sections; l++) {                            
        //                 #pragma omp atomic
        //                 GradSigma_3_tri_18N[k][i][j] += GradSigma_sec_tri_18_3[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
        //             }
        //         }
        //     }
        // }
        // // timing(time_s);
        // for (int i = 0; i < num_sections; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18; k++) {
        //             delete [] GradSigma_sec_tri_18_3[i][j][k];
        //         }
        //         delete [] GradSigma_sec_tri_18_3[i][j];
        //     }
        //     delete [] GradSigma_sec_tri_18_3[i];
        // }
        // delete [] GradSigma_sec_tri_18_3;  

        // // DivGradSigma
        // DivGradSigma_nu_sec_tri_18 = integral_divGradSigma_nu_dudv(G, num_sections, num_modes, nfp_ft, GLs, (int) GL_order, FuncDivGradSigma, \
        // R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, GradSigma_3_tri_18N, triangles_uv);
        // timing(time_s);
        // for (int i = 0; i < 3; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         delete [] GradSigma_3_tri_18N[i][j];
        //     }
        //     delete [] GradSigma_3_tri_18N[i];
        // }
        // delete [] GradSigma_3_tri_18N;
        // DivGradSigma_sec_tri_18 = solve_surfaces_sparse(DivGradSigma_nu_sec_tri_18, points_sequence, num_triangles, num_sections);
        // // timing(time_s);
        // for (int i = 0; i < num_sections; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         delete [] DivGradSigma_nu_sec_tri_18[i][j];
        //     }
        //     delete [] DivGradSigma_nu_sec_tri_18[i];
        // }
        // delete [] DivGradSigma_nu_sec_tri_18;  
        // // FT
        // DivGradSigma_tri_18N = new double *[num_triangles];
        // for (int i = 0; i < num_triangles; i++) {
        //     DivGradSigma_tri_18N[i] = new double [18*num_modes]();
        // }
        // #pragma omp parallel for
        // for (int i = 0; i < num_triangles; i++) {
        //     for (int j = 0; j < 18*num_modes; j++) {
        //         for (int l = 0; l < num_sections; l++) {                            
        //             #pragma omp atomic
        //             DivGradSigma_tri_18N[i][j] += DivGradSigma_sec_tri_18[l][i][j/num_modes] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
        //         }
        //     }
        // }
        // // timing(time_s);
        // for (int i = 0; i < num_sections; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         delete [] DivGradSigma_sec_tri_18[i][j];
        //     }
        //     delete [] DivGradSigma_sec_tri_18[i];
        // }
        // delete [] DivGradSigma_sec_tri_18;

        // // B2 for sigma
        // B2_DivGradSigma = Bpoisson_tri_18Nn(BG, DivGradSigma_tri_18N, num_triangles, num_modes, nfp_ft, error);
        // timing(time_s);
        // for (int i = 0; i < num_triangles; i++) {
        //     delete [] DivGradSigma_tri_18N[i];
        // }
        // delete [] DivGradSigma_tri_18N;

        // // solve poisson equation for sigma
        // double ***sigma_1_tri_18N = solve_pressure_sparse_LM(M2_poisson, B2_DivGradSigma, points_sequence, boundary_zeros, num_modes, num_triangles, error, time_s);
        // for (int i = 0; i < num_triangles; i++) {
        //     delete [] B2_DivGradSigma[i];
        // }
        // delete [] B2_DivGradSigma;
        // // print_p(sigma_1_tri_18N[0], G, triangles_uv, num_modes, nfp_ft);
        // // exit(0);

        // // Jperp
        // Jparallel_nu_sec_tri_18_3 = integral_Jparallel_nu_dudv(G, num_sections, num_modes, nfp_ft, GLs, GL_order, FuncJparallel, \
        //                             mf1, sigma_1_tri_18N[0], triangles_uv);
        // timing(time_s);
        // Jparallel_sec_tri_18_3 = solve_surfaces_sparse_num(Jparallel_nu_sec_tri_18_3, temp_sec_tri_18, B_sec_tri_18, points_sequence, num_sections, num_triangles, (int) 3);
        // // timing(time_s);
        // for (int i = 0; i < num_sections; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18; k++) {
        //             delete [] Jparallel_nu_sec_tri_18_3[i][j][k];
        //         }
        //         delete [] Jparallel_nu_sec_tri_18_3[i][j];
        //     }
        //     delete [] Jparallel_nu_sec_tri_18_3[i];
        // }
        // delete [] Jparallel_nu_sec_tri_18_3;  
        // // FT
        // Jparallel_3_tri_18N = new double **[3];
        // for (int i = 0; i < 3; i++) {
        //     Jparallel_3_tri_18N[i] = new double *[num_triangles];
        //     for (int j = 0; j < num_triangles; j++) {
        //         Jparallel_3_tri_18N[i][j] = new double [18*num_modes]();
        //     }
        // }
        // #pragma omp parallel for
        // for (int i = 0; i < num_triangles; i++) {
        //     for (int j = 0; j < 18*num_modes; j++) {
        //         for (int k = 0; k < 3; k++) {
        //             for (int l = 0; l < num_sections; l++) {                            
        //                 #pragma omp atomic
        //                 Jparallel_3_tri_18N[k][i][j] += Jparallel_sec_tri_18_3[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
        //             }
        //         }
        //     }
        // }
        // // timing(time_s);   
        // for (int i = 0; i < num_sections; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18; k++) {
        //             delete [] Jparallel_sec_tri_18_3[i][j][k];
        //         }
        //         delete [] Jparallel_sec_tri_18_3[i][j];
        //     }
        //     delete [] Jparallel_sec_tri_18_3[i];
        // }
        // delete [] Jparallel_sec_tri_18_3;  

        // //  J = Jperp + Jpara
        // J_3_tri_18N = new double **[3];
        // for (int i = 0; i < 3; i++) {
        //     J_3_tri_18N[i] = new double *[num_triangles];
        //     for (int j = 0; j < num_triangles; j++) {
        //         J_3_tri_18N[i][j] = new double [18*num_modes]();
        //     }
        // }
        // for (int i = 0; i < 3; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18*num_modes; k++) {
        //             // J_3_tri_18N[i][j][k] = Jperpendicular_3_tri_18N[i][j][k];
        //             J_3_tri_18N[i][j][k] = Jperpendicular_3_tri_18N[i][j][k] + Jparallel_3_tri_18N[i][j][k];
        //         }
        //     }
        // }
        // // print_B(Jperpendicular_3_tri_18N, G, triangles_uv, num_modes, nfp_ft, num_scalar);
        // // print_B(Jparallel_3_tri_18N, G, triangles_uv, num_modes, nfp_ft, num_scalar);
        // // print_B(J_3_tri_18N, G, triangles_uv, num_modes, nfp_ft, num_scalar);
        // // exit(0);
        // for (int i = 0; i < 3; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         delete [] Jperpendicular_3_tri_18N[i][j];
        //         delete [] Jparallel_3_tri_18N[i][j];
        //     }
        //     delete [] Jperpendicular_3_tri_18N[i];
        //     delete [] Jparallel_3_tri_18N[i];
        // }
        // delete [] Jperpendicular_3_tri_18N;
        // delete [] Jparallel_3_tri_18N;

        // // ---------------------------------------- iteration of B
        // // curl(J) * mu0
        // curlJ_nu_sec_tri_18_3 = integral_curlJ_nu_dudv(G, num_sections, num_modes, nfp_ft, GLs, (int) GL_order, FuncCurlJ, \
        // R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, J_3_tri_18N, triangles_uv);
        // timing(time_s); 
        // curlJ_sec_tri_18_3 = solve_surfaces_sparse_num(curlJ_nu_sec_tri_18_3, temp_sec_tri_18, B_sec_tri_18, points_sequence, num_sections, num_triangles, (int) 3);
        // // timing(time_s);
        // for (int i = 0; i < num_sections; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18; k++) {
        //             delete [] curlJ_nu_sec_tri_18_3[i][j][k];
        //         }
        //         delete [] curlJ_nu_sec_tri_18_3[i][j];
        //     }
        //     delete [] curlJ_nu_sec_tri_18_3[i];
        // }
        // delete [] curlJ_nu_sec_tri_18_3;  
        // // FT
        // curlJ_3_tri_18N = new double **[3];
        // for (int i = 0; i < 3; i++) {
        //     curlJ_3_tri_18N[i] = new double *[num_triangles];
        //     for (int j = 0; j < num_triangles; j++) {
        //         curlJ_3_tri_18N[i][j] = new double [18*num_modes]();
        //     }
        // }
        // #pragma omp parallel for
        // for (int i = 0; i < num_triangles; i++) {
        //     for (int j = 0; j < 18*num_modes; j++) {
        //         for (int k = 0; k < 3; k++) {
        //             for (int l = 0; l < num_sections; l++) {                            
        //                 #pragma omp atomic
        //                 curlJ_3_tri_18N[k][i][j] += curlJ_sec_tri_18_3[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
        //             }
        //         }
        //     }
        // }
        // // timing(time_s);   
        // for (int i = 0; i < num_sections; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18; k++) {
        //             delete [] curlJ_sec_tri_18_3[i][j][k];
        //         }
        //         delete [] curlJ_sec_tri_18_3[i][j];
        //     }
        //     delete [] curlJ_sec_tri_18_3[i];
        // }
        // delete [] curlJ_sec_tri_18_3;
        // // print_B(curlJ_3_tri_18N, G, triangles_uv, num_modes, nfp_ft, num_scalar);
        // // exit(0);

        // // B2 for curl(J)*mu0
        // B2_curlJ = new double *[num_triangles];
        // for (int i = 0; i < num_triangles; i++) {
        //     B2_curlJ[i] = new double [(18*num_modes)*num_scalar]();
        // }
        // #pragma omp parallel for
        // for (int i = 0; i < num_triangles; i++) {
        //     for (int k = 0; k < 18*num_modes; k++) {
        //         for (int t = 0; t < 18*num_modes; t++) {
        //             for (int p = 0; p < num_scalar; p++) {
        //                 #pragma omp atomic
        //                 B2_curlJ[i][k*num_scalar + p] += curlJ_3_tri_18N[p][i][t] * integral4(nfp, num_modes, t, k, 0, 0, error) * BG[i][t/num_modes][k/num_modes];
        //             }
        //         }   
        //     }
        // }
        // timing(time_s);
        // for (int i = 0; i < 3; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         delete [] curlJ_3_tri_18N[i][j];
        //     }
        //     delete [] curlJ_3_tri_18N[i];
        // }
        // delete [] curlJ_3_tri_18N;

        // // solve vector poisson equation for magnetic field
        // double ***mf2 = solve_magnetic_field_sparse_LM(M2, B2_curlJ, points_sequence, boundary_DoFs, num_modes, num_triangles, (int) num_scalar, error, time_s);
        // timing(time_s);
        // for (int i = 0; i < num_triangles; i++) {
        //     delete [] B2_curlJ[i];
        // }
        // delete [] B2_curlJ;

        // // print_B(mf2, G, triangles_uv, num_modes, nfp_ft, num_scalar); 

        // // mf.norm()
        // for (int i = 0; i < 3; i++) {
        //     norm1_B[i] = 0;
        //     norm2_B[i] = 0;
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18*num_modes; k++) {
        //             norm1_B[i] += mf1[i][j][k] * mf1[i][j][k];
        //             norm2_B[i] += mf2[i][j][k] * mf2[i][j][k];
        //         }
        //     }
        //     norm1_B[i] /= num_triangles * (18*num_modes);
        //     norm2_B[i] /= num_triangles * (18*num_modes);
        // }
        // printf("iter: %3d, mf, estimated error = %0.6e,  %0.6e,  %0.6e\n", iter + 1, abs(norm2_B[0]-norm1_B[0]) / norm1_B[0], abs(norm2_B[1]-norm1_B[1]) / norm1_B[1], abs(norm2_B[2]-norm1_B[2]) / norm1_B[2]);
        // fflush(stdout);

        // // exchange
        // for (int i = 0; i < 3; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18*num_modes; k++) {
        //             mf1[i][j][k] =  blending * mf2[i][j][k] + (1-blending) * mf1[i][j][k];
        //         }
        //     }
        // }
        // for (int i = 0; i < 3; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         delete [] mf2[i][j];
        //     }
        //     delete [] mf2[i];
        // }
        // delete [] mf2;

        // // upgrade for current density
        // #pragma omp parallel for collapse(2)
        // for (int i = 0; i < 3; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         for (int k = 0; k < 18*num_modes; k++) {
        //             J1[i][j][k] = J_3_tri_18N[i][j][k];
        //         }
        //     }
        // }
        // for (int i = 0; i < 3; i++) {
        //     for (int j = 0; j < num_triangles; j++) {
        //         delete [] J_3_tri_18N[i][j];
        //     }
        //     delete [] J_3_tri_18N[i];
        // }
        // delete [] J_3_tri_18N; 


    }
    printf("----------------------------------------------------\n\n\n\n\n\n");

    // recover pressure
    // div(grad(pressure)_parallel)
    divgradPparallel_nu_sec_tri_18 = integral_divJ_nu_dudv(G, num_sections, num_modes, nfp_ft, GLs, (int) GL_order, FuncDivGradPparallel, \
    R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, gradP1, triangles_uv);
    timing(time_s);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < num_triangles; j++) {
            delete [] gradP1[i][j];
        }
        delete [] gradP1[i];
    }
    delete [] gradP1;

    divgradPparallel_sec_tri_18 = solve_surfaces_sparse(divgradPparallel_nu_sec_tri_18, points_sequence, num_triangles, num_sections);


    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            delete [] divgradPparallel_nu_sec_tri_18[i][j];
        }
        delete [] divgradPparallel_nu_sec_tri_18[i];
    }
    delete [] divgradPparallel_nu_sec_tri_18;  
    // FT
    divgradPparallel_tri_18N = new double *[num_triangles];
    for (int i = 0; i < num_triangles; i++) {
        divgradPparallel_tri_18N[i] = new double [18*num_modes]();
    }
    #pragma omp parallel for
    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int l = 0; l < num_sections; l++) {                            
                #pragma omp atomic
                divgradPparallel_tri_18N[i][j] += divgradPparallel_sec_tri_18[l][i][j/num_modes] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) * dzeta;
            }
        }
    }
 
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles; j++) {
            delete [] divgradPparallel_sec_tri_18[i][j];
        }
        delete [] divgradPparallel_sec_tri_18[i];
    }
    delete [] divgradPparallel_sec_tri_18;
    // print_p(divgradPparallel_tri_18N, G, triangles_uv, num_modes, nfp_ft);

    // B2_poisson
    B2_divgradPperp = Bpoisson_tri_18Nn(BG, divgradPparallel_tri_18N, num_triangles, num_modes, nfp_ft, error);
    timing(time_s);
    for (int i = 0; i < num_triangles; i++) {
        delete [] divgradPparallel_tri_18N[i];
    }
    delete [] divgradPparallel_tri_18N;
    // solve poisson equation for pressure projection
    pres2 = solve_pressure_sparse_LM(M2_poisson, B2_divgradPperp, points_sequence, boundary_pressure, num_modes, num_triangles, error, time_s);
    // exchange
    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            pres1[i][j] = blending * pres2[0][i][j] + (1 - blending) * pres1[i][j];
        }
    }
    for (int i = 0; i < num_triangles; i++) {
        delete [] pres2[0][i];
    }
    delete [] pres2[0];
    delete [] pres2;

    print_p(pres1, G, triangles_uv, num_modes, nfp_ft);
    // print_B(mf1, G, triangles_uv, num_modes, nfp_ft, num_scalar);
    // print_B(J1, G, triangles_uv, num_modes, nfp_ft, num_scalar);
    
    

   
    printf("finished...\n");
    return 0;
}

