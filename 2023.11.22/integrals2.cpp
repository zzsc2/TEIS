#include "integrals.hpp"

//
double numerical_integration_pressure_over_triangle(double A, double B, double C, std::function<double(double , double , int , double **, vector<double>, int, double ***)> func, \
int position_tri, double **triangles_uv, vector<double> coefficients, int column, double ***G, \
double **GLs, int n) {
    // order = 5
    
    double c[25], x[25], y[25];
    for (int i = 0; i < 25; i++) {
        c[i] = GLs[0][i];
        x[i] = GLs[1][i];
        y[i] = GLs[2][i];
    }
    double uu, vv;
    double l1, l2, l3;
    double u, v;

    double h[25]={0}, val = 0;
    for (int k = 0; k < 25; k++) {

        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                h[k] += func(u, v, position_tri, triangles_uv, coefficients, column, G);
            }
        }

        for (int i = 0; i < 2 * n - 1; i++) {   
            for (int j = 0; j < 2 * n - i - 1; j++) {
                uu = (-x[k] + 2*(i-n) + 1) / (2*n);
                vv = (-y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                h[k] += func(u, v, position_tri, triangles_uv, coefficients, column, G);
            }
        }

        val += c[k] * h[k];  

    }

    val /= (4*n*n);
    val *= (A+B)*C/4;

    return val;
}

double ***integral_pressure_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double(double , double , int, double **, vector<double>)> functional, \
double **triangles_uv, vector<double> coefficients) {

    printf("calculating integral of pressure*nu_i...\n");
    int num_triangles_in_plane = (int) triangles_uv[0][6];
    
    std::function<double(double , double , int , double **, vector<double>, int, double ***)> functional_nu = \
    [&](double u, double v, int position_tri, double **triangles_uv, vector<double> coefficients, int column, double ***G) {

        double nu_i = 0, value;
        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }
        value = functional(u, v, position_tri, triangles_uv, coefficients) * nu_i;

        return value;
    };  

    double ***functional_nu_sec_tri_18 = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18[i] = new double *[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18[i][j] = new double [18];            
        }         
    }

    for (int i = 0; i < num_sections; i++) {
        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {
                // printf("i = %d, j = %d, k = %d\n", i, j, k);

                functional_nu_sec_tri_18[i][j][k] = numerical_integration_pressure_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                j, triangles_uv, coefficients, k, G,  \
                GLs, n);
            } 
        }
    }

    return functional_nu_sec_tri_18;
}


//
double ***numerical_integration_projection2_over_triangle(double A, double B, double C, std::function<double***(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int, int, int, double ***, double ***, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ***B_tri_3_18, double ***poisoon_2_tri_18, int column, \
double **GLs, int n) {
    // order = 5

    int num_derivates = 2;
    int num_scalar = 3;

    double c[25], x[25], y[25];
    for (int i = 0; i < 25; i++) {
        c[i] = GLs[0][i];
        x[i] = GLs[1][i];
        y[i] = GLs[2][i];
    }
    double uu, vv;
    double l1, l2, l3;
    double u, v;

    double ****h = new double ***[25];
    for (int i = 0; i < 25; i++) {
        h[i] = new double **[num_derivates];
        for (int j = 0; j < num_derivates; j++) {
            h[i][j] = new double *[num_scalar];
            for (int k = 0; k < num_scalar; k++) {
                h[i][j][k] = new double[num_scalar]();
            }
        }
    }

    double ***val = new double **[num_derivates];
    for (int i = 0; i < num_derivates; i++) {
        val[i] = new double *[num_scalar];
        for (int j = 0; j < num_scalar; j++) {
            val[i][j] = new double[num_scalar]();
        }
    }

    double ***ptr;
    for (int k = 0; k < 25; k++) {

        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri, num_modes, nfp_ft, num_sections, B_tri_3_18, poisoon_2_tri_18, column);
                for (int l = 0; l < num_derivates; l++) {
                    for (int s = 0 ; s < num_scalar; s++) {
                        for (int t = 0; t < num_scalar; t++) {
                            h[k][l][s][t] += ptr[l][s][t];
                        }
                    }
                }

                for (int l = 0; l < num_derivates; l++) {
                    for (int s = 0 ; s < num_scalar; s++) {
                        delete [] ptr[l][s];
                    }
                    delete [] ptr[l];
                }
                delete [] ptr;
                
            }
        }

        for (int i = 0; i < 2 * n - 1; i++) {   
            for (int j = 0; j < 2 * n - i - 1; j++) {
                uu = (-x[k] + 2*(i-n) + 1) / (2*n);
                vv = (-y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri, num_modes, nfp_ft, num_sections, B_tri_3_18, poisoon_2_tri_18, column);
                for (int l = 0; l < num_derivates; l++) {
                    for (int s = 0 ; s < num_scalar; s++) {
                        for (int t = 0; t < num_scalar; t++) {
                            h[k][l][s][t] += ptr[l][s][t];
                        }
                    }
                }

                for (int l = 0; l < num_derivates; l++) {
                    for (int s = 0 ; s < num_scalar; s++) {
                        delete [] ptr[l][s];
                    }
                    delete [] ptr[l];
                }
                delete [] ptr;
                
            }
        }

        for (int l = 0; l < num_derivates; l++) {
            for (int s = 0 ; s < num_scalar; s++) {
                for (int t = 0; t < num_scalar; t++) {
                    val[l][s][t] += c[k] * h[k][l][s][t];  
                }
            }
            
        }      

    }
    
    for (int i = 0; i < num_derivates; i++) {
        for (int s = 0 ; s < num_scalar; s++) {
            for (int t = 0; t < num_scalar; t++) {
                val[i][s][t] /= (4*n*n);
                val[i][s][t] *= (A+B)*C/4;
            }
        }
    }

    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < num_derivates; j++) {
            for (int k = 0; k < num_scalar; k++) {
                delete [] h[i][j][k];
            }
            delete [] h[i][j];
        }
        delete [] h[i];
    }
    delete [] h;

    return val;
}

double ******integral_projection2_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double ***(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int, int, int, double ***, double ***)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double **triangles_uv, int num_modes, int nfp_ft, double ***B_3_tri_18N, double ***poisoon_1_tri_18N) {

    printf("calculating integral of (Projection2)*nu_i...\n");
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double***(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int, int, int, double ***, double ***, int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ***B_tri_3_18, double ***poisoon_2_tri_18, int column) {

        double nu_i = 0;

        int num_scalar = 3;

        double ***value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri, num_modes, nfp_ft, num_sections, B_tri_3_18, poisoon_2_tri_18);
        for (int j = 0; j < num_scalar; j++) {
            for (int k = 0; k < num_scalar; k++) {
                value[0][j][k] *= nu_i;
                value[1][j][k] *= nu_i;
            }
        }

        return value;
    };

    double ******functional_nu_sec_tri_18_2_3_3 = new double *****[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_2_3_3[i] = new double ****[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_2_3_3[i][j] = new double ***[18];

        }
    }

    double ***B_tri_3_18 = new double **[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        B_tri_3_18[i] = new double *[3];
        for (int j = 0; j < 3; j++) {
            B_tri_3_18[i][j] = new double [18]();
        }
    }

    double ***poisoon_2_tri_18 = new double **[2];
    for (int i = 0; i < 2; i++) {
        poisoon_2_tri_18[i] = new double *[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            poisoon_2_tri_18[i][j] = new double [18]();
        }
    }

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18*num_modes; k++) {
                    #pragma omp atomic
                    B_tri_3_18[j][t][k/num_modes] += B_3_tri_18N[t][j][k] * fourier_basis_test(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
                }
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18*num_modes; k++) {
                #pragma omp atomic
                poisoon_2_tri_18[0][j][k/num_modes] += poisoon_1_tri_18N[0][j][k] *  fourier_basis_test(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
                #pragma omp atomic
                poisoon_2_tri_18[1][j][k/num_modes] += poisoon_1_tri_18N[0][j][k] * fourier_basis_test2(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {
                // printf("%d %d %d\n",i,j,k);
                // fflush(stdout);

                functional_nu_sec_tri_18_2_3_3[i][j][k] = numerical_integration_projection2_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, G, i, j, num_modes, nfp_ft, num_sections, B_tri_3_18, poisoon_2_tri_18, k, \
                                                                                             GLs, n);
            } 
        }

        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18; k++) {
                    B_tri_3_18[j][t][k] = 0.0;
                }
            }
            for (int k = 0; k < 18; k++) {
                poisoon_2_tri_18[0][j][k] = 0.0;
                poisoon_2_tri_18[1][j][k] = 0.0;
            }
        }

    }


    for (int j = 0; j < num_triangles_in_plane; j++) {
        for (int k = 0; k < 3; k++) {
            delete [] B_tri_3_18[j][k] ;
        }
        delete [] B_tri_3_18[j];
    }
    delete [] B_tri_3_18;

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            delete [] poisoon_2_tri_18[i][j];
        }
        delete [] poisoon_2_tri_18[i];
    }
    delete [] poisoon_2_tri_18;

    return functional_nu_sec_tri_18_2_3_3;
}





//
double *numerical_integration_gradP_over_triangle(double A, double B, double C, std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ***, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***pres_2_tri_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n) {
    // order = 5

    int num = 3;

    double c[25], x[25], y[25];
    for (int i = 0; i < 25; i++) {
        c[i] = GLs[0][i];
        x[i] = GLs[1][i];
        y[i] = GLs[2][i];
    }
    double uu, vv;
    double l1, l2, l3;
    double u, v;

    double h[25][num];
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < num; j++) {
            h[i][j] = 0;
        }
    }

    double *val = new double [num];
    for (int i = 0; i < num; i++) {
        val[i] = 0;
    }

    double *ptr;
    for (int k = 0; k < 25; k++) {

        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, pres_2_tri_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int i = 0; i < 2 * n - 1; i++) {   
            for (int j = 0; j < 2 * n - i - 1; j++) {
                uu = (-x[k] + 2*(i-n) + 1) / (2*n);
                vv = (-y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, pres_2_tri_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int l = 0; l < num; l++) {
            val[l] += c[k] * h[k][l];  
        }      

    }
    
    for (int i = 0; i < num; i++) {
        val[i] /= (4*n*n);
        val[i] *= (A+B)*C/4;
    }

    return val;
}

double ****integral_gradP_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double **pres_tri_18N, double **triangles_uv) {
    printf("calculating integral of (SMIE_grad(P))*nu_i...\n");
    fflush(stdout);
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ***, int , int , int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***pres_2_tri_18, int position_sec, int position_tri, \
    int column) {

        int num_derivaties = 3;

        double nu_i = 0;

        double *value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, G, pres_2_tri_18, position_sec, position_tri);
        for (int i = 0; i < num_derivaties; i++) {
            value[i] *= nu_i;
        }

        return value;
    };

    double ****functional_nu_sec_tri_18_3 = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_3[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_3[i][j] = new double *[18];
        }
    }

    double ***pres_2_tri_18 = new double **[2];
    for (int i = 0; i < 2; i++) {
        pres_2_tri_18[i] = new double *[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            pres_2_tri_18[i][j] = new double [18]();
        }
    }

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18*num_modes; k++) {
                #pragma omp atomic
                pres_2_tri_18[0][j][k/num_modes] += pres_tri_18N[j][k] *  fourier_basis_test(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
                #pragma omp atomic
                pres_2_tri_18[1][j][k/num_modes] += pres_tri_18N[j][k] * fourier_basis_test2(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_3[i][j][k] = numerical_integration_gradP_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, G, pres_2_tri_18, i, j, \
                                                                                             k, \
                                                                                             GLs, n);
            } 
        }

        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {
                pres_2_tri_18[0][j][k] = 0.0;
                pres_2_tri_18[1][j][k] = 0.0;
            }
        }

    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            delete [] pres_2_tri_18[i][j];
        }
        delete [] pres_2_tri_18[i];
    }
    delete [] pres_2_tri_18;

    return functional_nu_sec_tri_18_3;
}

//
double *numerical_integration_gradPparallel_over_triangle(double A, double B, double C, std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ***, double ***, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***mf_tri_3_18, double ***gradP_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n) {
    // order = 5

    int num = 3;

    double c[25], x[25], y[25];
    for (int i = 0; i < 25; i++) {
        c[i] = GLs[0][i];
        x[i] = GLs[1][i];
        y[i] = GLs[2][i];
    }
    double uu, vv;
    double l1, l2, l3;
    double u, v;

    double h[25][num];
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < num; j++) {
            h[i][j] = 0;
        }
    }

    double *val = new double [num];
    for (int i = 0; i < num; i++) {
        val[i] = 0;
    }

    double *ptr;
    for (int k = 0; k < 25; k++) {

        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, mf_tri_3_18, gradP_tri_3_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int i = 0; i < 2 * n - 1; i++) {   
            for (int j = 0; j < 2 * n - i - 1; j++) {
                uu = (-x[k] + 2*(i-n) + 1) / (2*n);
                vv = (-y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, mf_tri_3_18, gradP_tri_3_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int l = 0; l < num; l++) {
            val[l] += c[k] * h[k][l];  
        }      

    }
    
    for (int i = 0; i < num; i++) {
        val[i] /= (4*n*n);
        val[i] *= (A+B)*C/4;
    }

    return val;
}

double ****integral_gradPparallel_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***mf_3_tri_18N, double ***gradP_3_tri_18N, double **triangles_uv) {
    printf("calculating integral of (SMIE_grad(P) parallel)*nu_i...\n");
    fflush(stdout);
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ***, double ***, int , int , int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***mf_tri_3_18, double ***gradP_tri_3_18, int position_sec, int position_tri, \
    int column) {

        int num_derivaties = 3;

        double nu_i = 0;

        double *value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, G, mf_tri_3_18, gradP_tri_3_18, position_sec, position_tri);
        for (int i = 0; i < num_derivaties; i++) {
            value[i] *= nu_i;
        }

        return value;
    };

    double ****functional_nu_sec_tri_18_3 = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_3[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_3[i][j] = new double *[18];
        }
    }

    double ***mf_tri_3_18 = new double **[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        mf_tri_3_18[i] = new double *[3];
        for (int j = 0; j < 3; j++) {
            mf_tri_3_18[i][j] = new double [18]();
        }
    }

    double ***gradP_tri_3_18 = new double **[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        gradP_tri_3_18[i] = new double *[3];
        for (int j = 0; j < 3; j++) {
            gradP_tri_3_18[i][j] = new double [18]();
        }
    }

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18*num_modes; k++) {
                    #pragma omp atomic
                    mf_tri_3_18[j][t][k/num_modes] +=    mf_3_tri_18N[t][j][k] *    fourier_basis_test(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
                    #pragma omp atomic
                    gradP_tri_3_18[j][t][k/num_modes] += gradP_3_tri_18N[t][j][k] * fourier_basis_test(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
                }
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_3[i][j][k] = numerical_integration_gradPparallel_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, G, mf_tri_3_18, gradP_tri_3_18, i, j, \
                                                                                             k, \
                                                                                             GLs, n);
            } 
        }

        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18; k++) {
                    mf_tri_3_18[j][t][k] = 0.0;
                    gradP_tri_3_18[j][t][k] = 0.0;
                }
            }
        }

    }

    for (int j = 0; j < num_triangles_in_plane; j++) {
        for (int k = 0; k < 3; k++) {
            delete [] mf_tri_3_18[j][k];
            delete [] gradP_tri_3_18[j][k];
        }
        delete [] mf_tri_3_18[j];
        delete [] gradP_tri_3_18[j];
    }
    delete [] mf_tri_3_18;
    delete [] gradP_tri_3_18;

    return functional_nu_sec_tri_18_3;
}




//
double *numerical_integration_Jperpendicular_over_triangle(double A, double B, double C, std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ***, double ***, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***mf_tri_3_18, double ***gradP_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n) {
    // order = 5

    int num = 3;

    double c[25], x[25], y[25];
    for (int i = 0; i < 25; i++) {
        c[i] = GLs[0][i];
        x[i] = GLs[1][i];
        y[i] = GLs[2][i];
    }
    double uu, vv;
    double l1, l2, l3;
    double u, v;

    double h[25][num];
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < num; j++) {
            h[i][j] = 0;
        }
    }

    double *val = new double [num];
    for (int i = 0; i < num; i++) {
        val[i] = 0;
    }

    double *ptr;
    for (int k = 0; k < 25; k++) {

        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, mf_tri_3_18, gradP_tri_3_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int i = 0; i < 2 * n - 1; i++) {   
            for (int j = 0; j < 2 * n - i - 1; j++) {
                uu = (-x[k] + 2*(i-n) + 1) / (2*n);
                vv = (-y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, mf_tri_3_18, gradP_tri_3_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int l = 0; l < num; l++) {
            val[l] += c[k] * h[k][l];  
        }      

    }
    
    for (int i = 0; i < num; i++) {
        val[i] /= (4*n*n);
        val[i] *= (A+B)*C/4;
    }

    return val;
}

double ****integral_Jperpendicular_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***mf_3_tri_18N, double ***gradP_3_tri_18N, double **triangles_uv) {
    printf("calculating integral of (SMIE_X_perpendicular)*nu_i...\n");
    fflush(stdout);
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ***, double ***, int , int , int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***mf_tri_3_18, double ***gradP_tri_3_18, int position_sec, int position_tri, \
    int column) {

        int num_derivaties = 3;

        double nu_i = 0;

        double *value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, G, mf_tri_3_18, gradP_tri_3_18, position_sec, position_tri);
        for (int i = 0; i < num_derivaties; i++) {
            value[i] *= nu_i;
        }

        return value;
    };

    double ****functional_nu_sec_tri_18_3 = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_3[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_3[i][j] = new double *[18];
        }
    }

    double ***mf_tri_3_18 = new double **[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        mf_tri_3_18[i] = new double *[3];
        for (int j = 0; j < 3; j++) {
            mf_tri_3_18[i][j] = new double [18]();
        }
    }

    double ***gradP_tri_3_18 = new double **[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        gradP_tri_3_18[i] = new double *[3];
        for (int j = 0; j < 3; j++) {
            gradP_tri_3_18[i][j] = new double [18]();
        }
    }

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18*num_modes; k++) {
                    #pragma omp atomic
                    mf_tri_3_18[j][t][k/num_modes] += mf_3_tri_18N[t][j][k] * fourier_basis_test(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
                    #pragma omp atomic
                    gradP_tri_3_18[j][t][k/num_modes] += gradP_3_tri_18N[t][j][k] * fourier_basis_test(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
                }
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_3[i][j][k] = numerical_integration_Jperpendicular_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, G, mf_tri_3_18, gradP_tri_3_18, i, j, \
                                                                                             k, \
                                                                                             GLs, n);
            } 
        }

        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18; k++) {
                    mf_tri_3_18[j][t][k] = 0.0;
                    gradP_tri_3_18[j][t][k] = 0.0;
                }
            }
        }

    }

    for (int j = 0; j < num_triangles_in_plane; j++) {
        for (int k = 0; k < 3; k++) {
            delete [] mf_tri_3_18[j][k];
            delete [] gradP_tri_3_18[j][k];
        }
        delete [] mf_tri_3_18[j];
        delete [] gradP_tri_3_18[j];
    }
    delete [] mf_tri_3_18;
    delete [] gradP_tri_3_18;

    return functional_nu_sec_tri_18_3;
}

//
double *numerical_integration_curlJ_over_triangle(double A, double B, double C, std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, double ****, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****J_2_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n) {
    // order = 5

    int num = 3;

    double c[25], x[25], y[25];
    for (int i = 0; i < 25; i++) {
        c[i] = GLs[0][i];
        x[i] = GLs[1][i];
        y[i] = GLs[2][i];
    }
    double uu, vv;
    double l1, l2, l3;
    double u, v;

    double h[25][num];
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < num; j++) {
            h[i][j] = 0;
        }
    }

    double *val = new double [num];
    for (int i = 0; i < num; i++) {
        val[i] = 0;
    }

    double *ptr;
    for (int k = 0; k < 25; k++) {

        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int i = 0; i < 2 * n - 1; i++) {   
            for (int j = 0; j < 2 * n - i - 1; j++) {
                uu = (-x[k] + 2*(i-n) + 1) / (2*n);
                vv = (-y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int l = 0; l < num; l++) {
            val[l] += c[k] * h[k][l];  
        }      

    }
    
    for (int i = 0; i < num; i++) {
        val[i] /= (4*n*n);
        val[i] *= (A+B)*C/4;
    }

    return val;
}

double ****integral_curlJ_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ****, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***J_3_tri_18N, double **triangles_uv) {
    printf("calculating integral of (SMIE_curlJ)*nu_i...\n");
    fflush(stdout);
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ****, int , int , int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****J_2_tri_3_18, int position_sec, int position_tri, \
    int column) {

        int num_derivaties = 3;

        double nu_i = 0;

        double *value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, position_sec, position_tri);
        for (int i = 0; i < num_derivaties; i++) {
            value[i] *= nu_i;
        }

        return value;
    };

    double ****functional_nu_sec_tri_18_3 = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_3[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_3[i][j] = new double *[18];
        }
    }

    double ****J_2_tri_3_18 = new double ***[2];
    for (int i = 0; i < 2; i++) {
        J_2_tri_3_18[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            J_2_tri_3_18[i][j] = new double *[3];
            for (int k = 0; k < 3; k++) {
                J_2_tri_3_18[i][j][k] = new double [18]();

            }
        }
    }

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 18*num_modes; l++) {
                    #pragma omp atomic
                    J_2_tri_3_18[0][j][k][l/num_modes] += J_3_tri_18N[k][j][l] *  fourier_basis_test(l%num_modes, nfp_ft, (2*pi/num_sections) * i);
                    #pragma omp atomic
                    J_2_tri_3_18[1][j][k][l/num_modes] += J_3_tri_18N[k][j][l] * fourier_basis_test2(l%num_modes, nfp_ft, (2*pi/num_sections) * i);
                }
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_3[i][j][k] = numerical_integration_curlJ_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, i, j, \
                                                                                             k, \
                                                                                             GLs, n);
            } 
        }

        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18; k++) {
                    J_2_tri_3_18[0][j][t][k] = 0.0;
                    J_2_tri_3_18[1][j][t][k] = 0.0;
                }
            }
        }

    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                delete [] J_2_tri_3_18[i][j][k];
            }
            delete [] J_2_tri_3_18[i][j];
        }
        delete [] J_2_tri_3_18[i];
    }
    delete [] J_2_tri_3_18;

    return functional_nu_sec_tri_18_3;
}

//
double *numerical_integration_gradSigma_over_triangle(double A, double B, double C, std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, double ****, double ***, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****J_2_tri_3_18, double ***B_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n) {
    // order = 5

    int num = 3;

    double c[25], x[25], y[25];
    for (int i = 0; i < 25; i++) {
        c[i] = GLs[0][i];
        x[i] = GLs[1][i];
        y[i] = GLs[2][i];
    }
    double uu, vv;
    double l1, l2, l3;
    double u, v;

    double h[25][num];
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < num; j++) {
            h[i][j] = 0;
        }
    }

    double *val = new double [num];
    for (int i = 0; i < num; i++) {
        val[i] = 0;
    }

    double *ptr;
    for (int k = 0; k < 25; k++) {

        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, B_tri_3_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int i = 0; i < 2 * n - 1; i++) {   
            for (int j = 0; j < 2 * n - i - 1; j++) {
                uu = (-x[k] + 2*(i-n) + 1) / (2*n);
                vv = (-y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, B_tri_3_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int l = 0; l < num; l++) {
            val[l] += c[k] * h[k][l];  
        }      

    }
    
    for (int i = 0; i < num; i++) {
        val[i] /= (4*n*n);
        val[i] *= (A+B)*C/4;
    }

    return val;
}

double ****integral_gradSigma_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ****, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***J_3_tri_18N, double ***B_3_tri_18N, double **triangles_uv) {
    printf("calculating integral of (SMIE_gradSigma)*nu_i...\n");
    fflush(stdout);
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double*(double , double , double ***, double ***, double ***, double ***, double ***, double ****, double ***, int , int , int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****J_2_tri_3_18, double ***B_tri_3_18, int position_sec, int position_tri, \
    int column) {

        int num_derivaties = 3;

        double nu_i = 0;

        double *value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, B_tri_3_18, position_sec, position_tri);
        for (int i = 0; i < num_derivaties; i++) {
            value[i] *= nu_i;
        }

        return value;
    };

    double ****functional_nu_sec_tri_18_3 = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_3[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_3[i][j] = new double *[18];
        }
    }

    double ****J_2_tri_3_18 = new double ***[2];
    for (int i = 0; i < 2; i++) {
        J_2_tri_3_18[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            J_2_tri_3_18[i][j] = new double *[3];
            for (int k = 0; k < 3; k++) {
                J_2_tri_3_18[i][j][k] = new double [18]();

            }
        }
    }

    double ***B_tri_3_18 = new double **[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        B_tri_3_18[i] = new double *[3];
        for (int j = 0; j < 3; j++) {
            B_tri_3_18[i][j] = new double [18]();
        }
    }

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 18*num_modes; l++) {
                    #pragma omp atomic
                    J_2_tri_3_18[0][j][k][l/num_modes] += J_3_tri_18N[k][j][l] *  fourier_basis_test(l%num_modes, nfp_ft, (2*pi/num_sections) * i);
                    #pragma omp atomic
                    J_2_tri_3_18[1][j][k][l/num_modes] += J_3_tri_18N[k][j][l] * fourier_basis_test2(l%num_modes, nfp_ft, (2*pi/num_sections) * i);
                }
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18*num_modes; k++) {
                    #pragma omp atomic
                    B_tri_3_18[j][t][k/num_modes] += B_3_tri_18N[t][j][k] * fourier_basis_test(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
                }
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_3[i][j][k] = numerical_integration_gradSigma_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, B_tri_3_18, i, j, \
                                                                                             k, \
                                                                                             GLs, n);
            } 
        }

        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18; k++) {
                    J_2_tri_3_18[0][j][t][k] = 0.0;
                    J_2_tri_3_18[1][j][t][k] = 0.0;
                    B_tri_3_18[j][t][k] = 0.0;
                }
            }
        }

    }

    for (int j = 0; j < num_triangles_in_plane; j++) {
        for (int k = 0; k < 3; k++) {
            delete [] B_tri_3_18[j][k] ;
        }
        delete [] B_tri_3_18[j];
    }
    delete [] B_tri_3_18;

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                delete [] J_2_tri_3_18[i][j][k];
            }
            delete [] J_2_tri_3_18[i][j];
        }
        delete [] J_2_tri_3_18[i];
    }
    delete [] J_2_tri_3_18;

    return functional_nu_sec_tri_18_3;
}

//
double numerical_integration_divGradSigma_over_triangle(double A, double B, double C, std::function<double(double, double, double ***, double ***, double ***, double ***, double ***, double ****, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****GradSigma_2_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n) {
    
    // order = 5
    
    double c[25], x[25], y[25];
    for (int i = 0; i < 25; i++) {
        c[i] = GLs[0][i];
        x[i] = GLs[1][i];
        y[i] = GLs[2][i];
    }
    double uu, vv;
    double l1, l2, l3;
    double u, v;

    double h[25]={0}, val = 0;
    for (int k = 0; k < 25; k++) {

        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                h[k] += func(u, v, R, Z, Rzeta, Zzeta, G, GradSigma_2_tri_3_18, position_sec, position_tri, column);
            }
        }

        for (int i = 0; i < 2 * n - 1; i++) {   
            for (int j = 0; j < 2 * n - i - 1; j++) {
                uu = (-x[k] + 2*(i-n) + 1) / (2*n);
                vv = (-y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                h[k] += func(u, v, R, Z, Rzeta, Zzeta, G, GradSigma_2_tri_3_18, position_sec, position_tri, column);
            }
        }

        val += c[k] * h[k];  

    }

    val /= (4*n*n);
    val *= (A+B)*C/4;

    return val;
}

double ***integral_divGradSigma_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double (double, double, double ***, double ***, double ***, double ***, double ***, double ****, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***GradSigma_3_tri_18N, double **triangles_uv) {
    printf("calculating integral of (SMIE_divGradSigma)*nu_i...\n");
    fflush(stdout);
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double(double , double , double ***, double ***, double ***, double ***, double ***, double ****, int , int , int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****GradSigma_2_tri_3_18, int position_sec, int position_tri, \
    int column) {

        double nu_i = 0;

        double value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, G, GradSigma_2_tri_3_18, position_sec, position_tri) * nu_i;

        return value;
    };

    double ***functional_nu_sec_tri_18 = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18[i] = new double *[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18[i][j] = new double [18]();
        }
    }

    double ****GradSigma_2_tri_3_18 = new double ***[2];
    for (int i = 0; i < 2; i++) {
        GradSigma_2_tri_3_18[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            GradSigma_2_tri_3_18[i][j] = new double *[3];
            for (int k = 0; k < 3; k++) {
                GradSigma_2_tri_3_18[i][j][k] = new double [18]();

            }
        }
    }

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 18*num_modes; l++) {
                    #pragma omp atomic
                    GradSigma_2_tri_3_18[0][j][k][l/num_modes] += GradSigma_3_tri_18N[k][j][l] *  fourier_basis_test(l%num_modes, nfp_ft, (2*pi/num_sections) * i);
                    #pragma omp atomic
                    GradSigma_2_tri_3_18[1][j][k][l/num_modes] += GradSigma_3_tri_18N[k][j][l] * fourier_basis_test2(l%num_modes, nfp_ft, (2*pi/num_sections) * i);
                }
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18[i][j][k] = numerical_integration_divGradSigma_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, G, GradSigma_2_tri_3_18, i, j, \
                                                                                             k, \
                                                                                             GLs, n);
            } 
        }

        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18; k++) {
                    GradSigma_2_tri_3_18[0][j][t][k] = 0.0;
                    GradSigma_2_tri_3_18[1][j][t][k] = 0.0;
                }
            }
        }

    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                delete [] GradSigma_2_tri_3_18[i][j][k];
            }
            delete [] GradSigma_2_tri_3_18[i][j];
        }
        delete [] GradSigma_2_tri_3_18[i];
    }
    delete [] GradSigma_2_tri_3_18;

    return functional_nu_sec_tri_18;
}

//
double *numerical_integration_Jparallel_over_triangle(double A, double B, double C, std::function<double*(double, double, double ***, double ***, double **, int , int , int)> func, \
double ***G, double ***mf_3_tri_18, double **sigma_tri_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n) {
    // order = 5

    int num = 3;

    double c[25], x[25], y[25];
    for (int i = 0; i < 25; i++) {
        c[i] = GLs[0][i];
        x[i] = GLs[1][i];
        y[i] = GLs[2][i];
    }
    double uu, vv;
    double l1, l2, l3;
    double u, v;

    double h[25][num];
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < num; j++) {
            h[i][j] = 0;
        }
    }

    double *val = new double [num];
    for (int i = 0; i < num; i++) {
        val[i] = 0;
    }

    double *ptr;
    for (int k = 0; k < 25; k++) {

        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, G, mf_3_tri_18, sigma_tri_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int i = 0; i < 2 * n - 1; i++) {   
            for (int j = 0; j < 2 * n - i - 1; j++) {
                uu = (-x[k] + 2*(i-n) + 1) / (2*n);
                vv = (-y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                ptr = func(u, v, G, mf_3_tri_18, sigma_tri_18, position_sec, position_tri, column);
                for (int l = 0; l < num; l++) {
                    h[k][l] += ptr[l];
                }
                delete [] ptr;
            }
        }

        for (int l = 0; l < num; l++) {
            val[l] += c[k] * h[k][l];  
        }      

    }
    
    for (int i = 0; i < num; i++) {
        val[i] /= (4*n*n);
        val[i] *= (A+B)*C/4;
    }

    return val;
}

double ****integral_Jparallel_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double **, int, int)> functional, \
double ***mf_3_tri_18N, double **sigma_tri_18N, double **triangles_uv) {
    printf("calculating integral of (SMIE_X_parallel)*nu_i...\n");
    fflush(stdout);
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double*(double , double , double ***, double ***, double **, int , int , int)> functional_nu = \
    [&](double u, double v, double ***G, double ***mf_3_tri_18, double **sigma_tri_18, int position_sec, int position_tri, \
    int column) {

        int num_derivaties = 3;

        double nu_i = 0;

        double *value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, G, mf_3_tri_18, sigma_tri_18, position_sec, position_tri);
        for (int i = 0; i < num_derivaties; i++) {
            value[i] *= nu_i;
        }

        return value;
    };

    double ****functional_nu_sec_tri_18_3 = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_3[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_3[i][j] = new double *[18];
        }
    }

    double ***mf_tri_3_18 = new double **[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        mf_tri_3_18[i] = new double *[3];
        for (int j = 0; j < 3; j++) {
            mf_tri_3_18[i][j] = new double [18]();
        }
    }

    double **sigma_tri_18 = new double *[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        sigma_tri_18[i] = new double [18]();
    }

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18*num_modes; k++) {
                    #pragma omp atomic
                    mf_tri_3_18[j][t][k/num_modes] += mf_3_tri_18N[t][j][k] * fourier_basis_test(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
                }
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18*num_modes; k++) {
                #pragma omp atomic
                sigma_tri_18[j][k/num_modes] += sigma_tri_18N[j][k] *  fourier_basis_test(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_3[i][j][k] = numerical_integration_Jparallel_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             G, mf_tri_3_18, sigma_tri_18, i, j, \
                                                                                             k, \
                                                                                             GLs, n);
            } 
        }

        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18; k++) {
                    mf_tri_3_18[j][t][k] = 0.0;
                }
            }
            for (int k = 0; k < 18; k++) {
                sigma_tri_18[j][k] = 0.0;
            }
        }

    }

    for (int j = 0; j < num_triangles_in_plane; j++) {
        for (int k = 0; k < 3; k++) {
            delete [] mf_tri_3_18[j][k] ;
        }
        delete [] mf_tri_3_18[j];
    }
    delete [] mf_tri_3_18;

    for (int i = 0; i < num_triangles_in_plane; i++) {
        delete [] sigma_tri_18[i];
    }
    delete [] sigma_tri_18;

    return functional_nu_sec_tri_18_3;
}






//
double numerical_integration_pres_over_triangle(double A, double B, double C, std::function<double(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, double ***, double ****, double ***, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, double ***pressure_3_tri_18, double ****GradPparallel_2_tri_3_18, double ***B_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n) {
    
    // order = 5
    
    double c[25], x[25], y[25];
    for (int i = 0; i < 25; i++) {
        c[i] = GLs[0][i];
        x[i] = GLs[1][i];
        y[i] = GLs[2][i];
    }
    double uu, vv;
    double l1, l2, l3;
    double u, v;

    double h[25]={0}, val = 0;
    for (int k = 0; k < 25; k++) {

        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                h[k] += func(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, pressure_3_tri_18, GradPparallel_2_tri_3_18, B_tri_3_18, position_sec, position_tri, column);
            }
        }

        for (int i = 0; i < 2 * n - 1; i++) {   
            for (int j = 0; j < 2 * n - i - 1; j++) {
                uu = (-x[k] + 2*(i-n) + 1) / (2*n);
                vv = (-y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                h[k] += func(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, pressure_3_tri_18, GradPparallel_2_tri_3_18, B_tri_3_18, position_sec, position_tri, column);
            }
        }

        val += c[k] * h[k];  

    }

    val /= (4*n*n);
    val *= (A+B)*C/4;

    return val;
}

double ***integral_pres_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double (double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, double ***, double ****, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double **pressure_tri_18N, double ***GradPparallel_tri_3_18N, double ***B_3_tri_18N, double **triangles_uv) {
    printf("calculating integral of (SMIE_pressure)*nu_i...\n");
    fflush(stdout);
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double(double , double , double ***, double ***, double ***, double ***, double ***, double ***, double ***, double ***, double ****, double ***, int , int , int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, double ***pressure_3_tri_18, double ****GradPparallel_2_tri_3_18, double ***B_tri_3_18, int position_sec, int position_tri, \
    int column) {

        double nu_i = 0;

        double value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, pressure_3_tri_18, GradPparallel_2_tri_3_18, B_tri_3_18, position_sec, position_tri) * nu_i;

        return value;
    };

    double ***functional_nu_sec_tri_18 = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18[i] = new double *[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18[i][j] = new double [18]();
        }
    }

    double ***pressure_3_tri_18 = new double **[3];
    for (int i = 0; i < 3; i++) {
        pressure_3_tri_18[i] = new double *[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            pressure_3_tri_18[i][j] = new double [18]();
        }
    }

    double ***B_tri_3_18 = new double **[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        B_tri_3_18[i] = new double *[3];
        for (int j = 0; j < 3; j++) {
            B_tri_3_18[i][j] = new double [18]();
        }
    }

    double ****GradPparallel_2_tri_3_18 = new double ***[2];
    for (int i = 0; i < 2; i++) {
        GradPparallel_2_tri_3_18[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            GradPparallel_2_tri_3_18[i][j] = new double *[3];
            for (int k = 0; k < 3; k++) {
                GradPparallel_2_tri_3_18[i][j][k] = new double [18]();
            }
        }
    }

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 18*num_modes; l++) {
                    #pragma omp atomic
                    GradPparallel_2_tri_3_18[0][j][k][l/num_modes] += GradPparallel_tri_3_18N[k][j][l] *  fourier_basis_test(l%num_modes, nfp_ft, (2.0*pi/num_sections) * i);
                    #pragma omp atomic
                    GradPparallel_2_tri_3_18[1][j][k][l/num_modes] += GradPparallel_tri_3_18N[k][j][l] * fourier_basis_test2(l%num_modes, nfp_ft, (2.0*pi/num_sections) * i);
                    #pragma omp atomic
                    B_tri_3_18[j][k][l/num_modes] += B_3_tri_18N[k][j][l] *  fourier_basis_test(l%num_modes, nfp_ft, (2.0*pi/num_sections) * i);
                }
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int l = 0; l < 18*num_modes; l++) {
                #pragma omp atomic
                pressure_3_tri_18[0][j][l/num_modes] += pressure_tri_18N[j][l] *  fourier_basis_test(l%num_modes, nfp_ft, (2.0*pi/num_sections) * i);
                #pragma omp atomic
                pressure_3_tri_18[1][j][l/num_modes] += pressure_tri_18N[j][l] * fourier_basis_test2(l%num_modes, nfp_ft, (2.0*pi/num_sections) * i);
                #pragma omp atomic
                pressure_3_tri_18[2][j][l/num_modes] += pressure_tri_18N[j][l] * fourier_basis_test3(l%num_modes, nfp_ft, (2.0*pi/num_sections) * i);
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18[i][j][k] = numerical_integration_pres_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, pressure_3_tri_18, GradPparallel_2_tri_3_18, B_tri_3_18, i, j, \
                                                                                             k, \
                                                                                             GLs, n);
            } 
        }

        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18; k++) {
                    GradPparallel_2_tri_3_18[0][j][t][k] = 0.0;
                    GradPparallel_2_tri_3_18[1][j][t][k] = 0.0;
                    B_tri_3_18[j][t][k] = 0.0;
                    pressure_3_tri_18[t][j][k] = 0.0;
                }
            }
        }

    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                delete [] GradPparallel_2_tri_3_18[i][j][k];
            }
            delete [] GradPparallel_2_tri_3_18[i][j];
        }
        delete [] GradPparallel_2_tri_3_18[i];
    }
    delete [] GradPparallel_2_tri_3_18;

    for (int j = 0; j < num_triangles_in_plane; j++) {
        for (int k = 0; k < 3; k++) {
            delete [] B_tri_3_18[j][k] ;
        }
        delete [] B_tri_3_18[j];
    }
    delete [] B_tri_3_18;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            delete [] pressure_3_tri_18[i][j];
        }
        delete [] pressure_3_tri_18[i];
    }
    delete [] pressure_3_tri_18;

    return functional_nu_sec_tri_18;
}











//
double numerical_integration_divJ_over_triangle(double A, double B, double C, std::function<double(double, double, double ***, double ***, double ***, double ***, double ***, double ****, int , int , int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****J_2_tri_3_18, int position_sec, int position_tri, \
int column, \
double **GLs, int n) {
    
    // order = 5
    
    double c[25], x[25], y[25];
    for (int i = 0; i < 25; i++) {
        c[i] = GLs[0][i];
        x[i] = GLs[1][i];
        y[i] = GLs[2][i];
    }
    double uu, vv;
    double l1, l2, l3;
    double u, v;

    double h[25]={0}, val = 0;
    for (int k = 0; k < 25; k++) {

        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                h[k] += func(u, v, R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, position_sec, position_tri, column);
            }
        }

        for (int i = 0; i < 2 * n - 1; i++) {   
            for (int j = 0; j < 2 * n - i - 1; j++) {
                uu = (-x[k] + 2*(i-n) + 1) / (2*n);
                vv = (-y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                h[k] += func(u, v, R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, position_sec, position_tri, column);
            }
        }

        val += c[k] * h[k];  

    }

    val /= (4*n*n);
    val *= (A+B)*C/4;

    return val;
}

double ***integral_divJ_nu_dudv(double ***G, int num_sections, int num_modes, int nfp_ft, double **GLs, int n, \
std::function<double (double, double, double ***, double ***, double ***, double ***, double ***, double ****, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***J_3_tri_18N, double **triangles_uv) {
    printf("calculating integral of (SMIE_divJ)*nu_i...\n");
    fflush(stdout);
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double(double , double , double ***, double ***, double ***, double ***, double ***, double ****, int , int , int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****J_2_tri_3_18, int position_sec, int position_tri, \
    int column) {

        int num_derivaties = 3;

        double nu_i = 0;

        double value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, position_sec, position_tri) * nu_i;

        return value;
    };

    double ***functional_nu_sec_tri_18 = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18[i] = new double *[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18[i][j] = new double [18]();
        }
    }

    double ****J_2_tri_3_18 = new double ***[2];
    for (int i = 0; i < 2; i++) {
        J_2_tri_3_18[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            J_2_tri_3_18[i][j] = new double *[3];
            for (int k = 0; k < 3; k++) {
                J_2_tri_3_18[i][j][k] = new double [18]();
            }
        }
    }

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 18*num_modes; l++) {
                    #pragma omp atomic
                    J_2_tri_3_18[0][j][k][l/num_modes] += J_3_tri_18N[k][j][l] *  fourier_basis_test(l%num_modes, nfp_ft, (2*pi/num_sections) * i);
                    #pragma omp atomic
                    J_2_tri_3_18[1][j][k][l/num_modes] += J_3_tri_18N[k][j][l] * fourier_basis_test2(l%num_modes, nfp_ft, (2*pi/num_sections) * i);
                }
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18[i][j][k] = numerical_integration_divJ_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, G, J_2_tri_3_18, i, j, \
                                                                                             k, \
                                                                                             GLs, n);
            } 
        }

        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18; k++) {
                    J_2_tri_3_18[0][j][t][k] = 0.0;
                    J_2_tri_3_18[1][j][t][k] = 0.0;
                }
            }
        }

    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                delete [] J_2_tri_3_18[i][j][k];
            }
            delete [] J_2_tri_3_18[i][j];
        }
        delete [] J_2_tri_3_18[i];
    }
    delete [] J_2_tri_3_18;

    return functional_nu_sec_tri_18;
}





