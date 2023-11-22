#include "meshes.hpp"

double ***triangles_mapping(double **triangles_xy, double **fourier_coefficients, int num_sections, int nfp) {
    // linear mapping
    printf("mapping...\n");
    
    int num_triangles = (int) triangles_xy[0][6];
    int num_coefficients = (int) fourier_coefficients[0][6];

    double **control_points_u = new double *[num_triangles];
    double **control_points_v = new double *[num_triangles];
    
    double **theta_cp_u = new double *[num_triangles];
    double **theta_cp_v = new double *[num_triangles];

    double **p = new double *[num_triangles];
    double **q = new double *[num_triangles];

    for (int i = 0; i < num_triangles; i++) {
        
        control_points_u[i] = new double[12];
        control_points_v[i] = new double[12];

        theta_cp_u[i] = new double[6];
        theta_cp_v[i] = new double[6];

        p[i] = new double[3];
        q[i] = new double[3];

        for (int j = 0; j < 12; j+=4) {

            control_points_u[i][j  ] =-sqrt(1 - triangles_xy[i][j/2+1] * triangles_xy[i][j/2+1]);
            control_points_u[i][j+1] = triangles_xy[i][j/2+1];
            control_points_u[i][j+2] = sqrt(1 - triangles_xy[i][j/2+1] * triangles_xy[i][j/2+1]);
            control_points_u[i][j+3] = triangles_xy[i][j/2+1];

            control_points_v[i][j  ] = triangles_xy[i][j/2];
            control_points_v[i][j+1] =-sqrt(1 - triangles_xy[i][j/2] * triangles_xy[i][j/2]);
            control_points_v[i][j+2] = triangles_xy[i][j/2];
            control_points_v[i][j+3] = sqrt(1 - triangles_xy[i][j/2] * triangles_xy[i][j/2]);

            theta_cp_u[i][j/2  ] = atan2(control_points_u[i][j+1], control_points_u[i][j+0]);
            theta_cp_u[i][j/2+1] = atan2(control_points_u[i][j+3], control_points_u[i][j+2]);
            theta_cp_v[i][j/2  ] = atan2(control_points_v[i][j+1], control_points_v[i][j+0]);  
            theta_cp_v[i][j/2+1] = atan2(control_points_v[i][j+3], control_points_v[i][j+2]);

            if (control_points_u[i][j+0] - control_points_u[i][j+2] != 0) {
                p[i][j/4] = (triangles_xy[i][j/2+0] - control_points_u[i][j+2]) / (control_points_u[i][j+0] - control_points_u[i][j+2]);
            }
            else {
                p[i][j/4] = 0;
            }

            if (control_points_v[i][j+1] - control_points_v[i][j+3] != 0)  {
                q[i][j/4] = (triangles_xy[i][j/2+1] - control_points_v[i][j+3]) / (control_points_v[i][j+1] - control_points_v[i][j+3]);
            }
            else {
                q[i][j/4] = 0;
            }
        }
    }
   
    double *zeta_nfp = new double [num_sections];
    double ***triangles_RZ = new double **[num_sections];
    // #pragma omp parallel for
    for (int i = 0; i < num_sections; i++) {

        triangles_RZ[i] = new double *[num_triangles];
        zeta_nfp[i] = (2*pi)/num_sections*i*nfp;

        for (int j = 0; j < num_triangles; j++) {

            triangles_RZ[i][j] = new double[6];

            triangles_RZ[i][j][0] = p[j][0] * boundary_R(theta_cp_u[j][0], zeta_nfp[i], fourier_coefficients) + (1 - p[j][0]) * boundary_R(theta_cp_u[j][1], zeta_nfp[i], fourier_coefficients);
            triangles_RZ[i][j][1] = q[j][0] * boundary_Z(theta_cp_v[j][0], zeta_nfp[i], fourier_coefficients) + (1 - q[j][0]) * boundary_Z(theta_cp_v[j][1], zeta_nfp[i], fourier_coefficients);

            triangles_RZ[i][j][2] = p[j][1] * boundary_R(theta_cp_u[j][2], zeta_nfp[i], fourier_coefficients) + (1 - p[j][1]) * boundary_R(theta_cp_u[j][3], zeta_nfp[i], fourier_coefficients);
            triangles_RZ[i][j][3] = q[j][1] * boundary_Z(theta_cp_v[j][2], zeta_nfp[i], fourier_coefficients) + (1 - q[j][1]) * boundary_Z(theta_cp_v[j][3], zeta_nfp[i], fourier_coefficients);

            triangles_RZ[i][j][4] = p[j][2] * boundary_R(theta_cp_u[j][4], zeta_nfp[i], fourier_coefficients) + (1 - p[j][2]) * boundary_R(theta_cp_u[j][5], zeta_nfp[i], fourier_coefficients);
            triangles_RZ[i][j][5] = q[j][2] * boundary_Z(theta_cp_v[j][4], zeta_nfp[i], fourier_coefficients) + (1 - q[j][2]) * boundary_Z(theta_cp_v[j][5], zeta_nfp[i], fourier_coefficients);

        }
    }

    triangles_RZ[0][0] = enlarge_1d(triangles_RZ[0][0], 6, 7);
    triangles_RZ[0][0][6] = (int) triangles_xy[0][6];  
    for (int i = 0; i < num_triangles; i++) {
        delete [] control_points_u[i];
        delete [] control_points_v[i]; 
        delete [] theta_cp_u[i]; 
        delete [] theta_cp_v[i];
        delete [] p[i];
        delete [] q[i]; 
    }      
    delete [] control_points_u;
    delete [] control_points_v;
    delete [] theta_cp_u;
    delete [] theta_cp_v;
    delete [] p;
    delete [] q;
    delete [] zeta_nfp;
    
    // int j = 20;
    // for (int i = 0; i < num_triangles; i++) {
    //     printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n",triangles_RZ[j][i][0],triangles_RZ[j][i][1],triangles_RZ[j][i][2],triangles_RZ[j][i][3],triangles_RZ[j][i][4],triangles_RZ[j][i][5]);
    // }
    // exit(0);

    // for (int i = 0; i < num_triangles; i++) {
    //     printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n",triangles_xy[i][0],triangles_xy[i][1],triangles_xy[i][2],triangles_xy[i][3],triangles_xy[i][4],triangles_xy[i][5]);
    // }
    // exit(0);

    return triangles_RZ;
}

double coordinate_R(double u, double v, double x0, double y0, double zeta, double **fourier_coefficients, int nfp, double sine, double cosine) {
    int num_coefficients = (int) fourier_coefficients[0][6];

    double x = x0 + u*cosine - v*sine;
    double y = y0 + u*sine + v*cosine;
    double xx;
    if (abs(y) > 1.0) {
        xx = 0;
    }
    else {
        xx = sqrt(1-y*y);
    }
    double theta1 = atan2(y, -xx);
    double theta3 = atan2(y,  xx);

    double p = 1.0;
    if (xx != 0.0) {
        p = x / (2*xx) + 0.5;
    }

    return p * boundary_R(theta3, zeta*nfp, fourier_coefficients) + (1-p) * boundary_R(theta1, zeta*nfp, fourier_coefficients);
}

double coordinate_Z(double u, double v, double x0, double y0, double zeta, double **fourier_coefficients, int nfp, double sine, double cosine) {
    int num_coefficients = (int) fourier_coefficients[0][6];

    double x = x0 + u*cosine - v*sine;
    double y = y0 + u*sine + v*cosine;
    double yy;
    if (abs(x) > 1.0) {
        yy = 0;
    }
    else {
        yy = sqrt(1-x*x);
    }
    double theta2 = atan2(-yy, x);
    double theta4 = atan2( yy, x);
    
    double q = 1.0;
    if (yy != 0.0) {
        q = y / (2*yy) + 0.5;
    }

    return q * boundary_Z(theta4, zeta*nfp, fourier_coefficients) + (1-q) * boundary_Z(theta2, zeta*nfp, fourier_coefficients);
}

double boundary_R(double theta, double zeta_nfp, double **fourier_coefficients) {
    int num_coefficients = (int) fourier_coefficients[0][6];
    
    double value = 0;
    double angle = 0;
    for (int i = 0; i < num_coefficients; i++) {  
        angle = (fourier_coefficients[i][1]*theta + fourier_coefficients[i][0]*zeta_nfp);
        value += fourier_coefficients[i][2] * cosine_lookup_approximation(angle) + fourier_coefficients[i][4] * sine_lookup_approximation(angle);
        // value += fourier_coefficients[i][2] * cos(angle) + fourier_coefficients[i][4] * sin(angle);
    }
    return value;
}

double boundary_Z(double theta, double zeta_nfp, double **fourier_coefficients) {
    int num_coefficients = (int) fourier_coefficients[0][6];

    double value = 0;
    double angle = 0;
    for (int i = 0; i < num_coefficients; i++) {
        angle = (fourier_coefficients[i][1]*theta + fourier_coefficients[i][0]*zeta_nfp);
        value += fourier_coefficients[i][3] * sine_lookup_approximation(angle) + fourier_coefficients[i][5] * cosine_lookup_approximation(angle); 
        // value += fourier_coefficients[i][3] * sin(angle) + fourier_coefficients[i][5] * cos(angle); 
    }

    return value;
}

void triangles_reorder(double ***triangles_RZ, double ** triangles_xy, int num_sections) {
    printf("rearranging triangles...\n");
    int num_triangles_in_plane = (int) triangles_RZ[0][0][6];

    double l_12, l_23, l_31;
    double dx1, dy1, dx2, dy2, dx3, dy3;
    double ex1[2], ex2[2], ex3[2];
    for (int j = 0; j < num_triangles_in_plane; j++) {
        dx1 = triangles_xy[j][0]-triangles_xy[j][2];
        dy1 = triangles_xy[j][1]-triangles_xy[j][3];
        dx2 = triangles_xy[j][2]-triangles_xy[j][4];
        dy2 = triangles_xy[j][3]-triangles_xy[j][5];
        dx3 = triangles_xy[j][4]-triangles_xy[j][0];
        dy3 = triangles_xy[j][5]-triangles_xy[j][1];

        l_12 = dx1 * dx1 + dy1 * dy1;
        l_23 = dx2 * dx2 + dy2 * dy2;
        l_31 = dx3 * dx3 + dy3 * dy3;

        if (l_23 > l_12) {
            ex1[0] = triangles_xy[j][0];
            ex1[1] = triangles_xy[j][1];
            ex2[0] = triangles_xy[j][2];
            ex2[1] = triangles_xy[j][3];
            ex3[0] = triangles_xy[j][4];
            ex3[1] = triangles_xy[j][5];
            if (l_23 > l_31) {
                triangles_xy[j][0] = ex2[0];
                triangles_xy[j][1] = ex2[1];
                triangles_xy[j][2] = ex3[0];
                triangles_xy[j][3] = ex3[1];
                triangles_xy[j][4] = ex1[0];
                triangles_xy[j][5] = ex1[1];
                for (int i = 0; i < num_sections; i++) {
                    ex1[0] = triangles_RZ[i][j][0];
                    ex1[1] = triangles_RZ[i][j][1];
                    ex2[0] = triangles_RZ[i][j][2];
                    ex2[1] = triangles_RZ[i][j][3];
                    ex3[0] = triangles_RZ[i][j][4];
                    ex3[1] = triangles_RZ[i][j][5];
                    triangles_RZ[i][j][0] = ex2[0];
                    triangles_RZ[i][j][1] = ex2[1];
                    triangles_RZ[i][j][2] = ex3[0];
                    triangles_RZ[i][j][3] = ex3[1];
                    triangles_RZ[i][j][4] = ex1[0];
                    triangles_RZ[i][j][5] = ex1[1];
                }
            }
            else {
                triangles_xy[j][0] = ex3[0];
                triangles_xy[j][1] = ex3[1];
                triangles_xy[j][2] = ex1[0];
                triangles_xy[j][3] = ex1[1];
                triangles_xy[j][4] = ex2[0];
                triangles_xy[j][5] = ex2[1];
                for (int i = 0; i < num_sections; i++) {
                    ex1[0] = triangles_RZ[i][j][0];
                    ex1[1] = triangles_RZ[i][j][1];
                    ex2[0] = triangles_RZ[i][j][2];
                    ex2[1] = triangles_RZ[i][j][3];
                    ex3[0] = triangles_RZ[i][j][4];
                    ex3[1] = triangles_RZ[i][j][5];
                    triangles_RZ[i][j][0] = ex3[0];
                    triangles_RZ[i][j][1] = ex3[1];
                    triangles_RZ[i][j][2] = ex1[0];
                    triangles_RZ[i][j][3] = ex1[1];
                    triangles_RZ[i][j][4] = ex2[0];
                    triangles_RZ[i][j][5] = ex2[1];
                }
            }
        }
        else if (l_12 < l_31) {            
            ex1[0] = triangles_xy[j][0];
            ex1[1] = triangles_xy[j][1];
            ex2[0] = triangles_xy[j][2];
            ex2[1] = triangles_xy[j][3];
            ex3[0] = triangles_xy[j][4];
            ex3[1] = triangles_xy[j][5];

            triangles_xy[j][0] = ex3[0];
            triangles_xy[j][1] = ex3[1];
            triangles_xy[j][2] = ex1[0];
            triangles_xy[j][3] = ex1[1];
            triangles_xy[j][4] = ex2[0];
            triangles_xy[j][5] = ex2[1];
            for (int i = 0; i < num_sections; i++) {
                ex1[0] = triangles_RZ[i][j][0];
                ex1[1] = triangles_RZ[i][j][1];
                ex2[0] = triangles_RZ[i][j][2];
                ex2[1] = triangles_RZ[i][j][3];
                ex3[0] = triangles_RZ[i][j][4];
                ex3[1] = triangles_RZ[i][j][5];
                triangles_RZ[i][j][0] = ex3[0];
                triangles_RZ[i][j][1] = ex3[1];
                triangles_RZ[i][j][2] = ex1[0];
                triangles_RZ[i][j][3] = ex1[1];
                triangles_RZ[i][j][4] = ex2[0];
                triangles_RZ[i][j][5] = ex2[1];
            }                
        }
    }

    // int pos = 0;
    // for (int i = 0; i < num_triangles_in_plane; i++) {
    //     printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n",triangles_RZ[pos][i][0],triangles_RZ[pos][i][1],triangles_RZ[pos][i][2],triangles_RZ[pos][i][3],triangles_RZ[pos][i][4],triangles_RZ[pos][i][5]);
    // }
    // printf("\n\n\n");
    // exit(0);
}

int **point_retrieval(double **triangles_xy, int num_edges, int num_vertices, double error) {
    printf("retrieving points...\n");

    int num_triangles = (int) triangles_xy[0][6];

    int num1 = 0;
    double **p_internal = new double *[num_vertices-num_edges];
    for (int i = 0; i < num_vertices-num_edges; i++) {
        p_internal[i] = new double [2];
        for (int j = 0; j < 2; j++) {
            p_internal[i][j] = 9999.0;
        }
    }

    int num2 = 0;
    double **p_boundary = new double *[num_edges];
    for (int i = 0; i < num_edges; i++) {
        p_boundary[i] = new double [2];
        for (int j = 0; j < 2; j++) {
            p_boundary[i][j] = 9999.0;
        }
    }

    int num3 = 0;
    double **p_global = new double *[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        p_global[i] = new double [2];
        for (int j = 0; j < 2; j++) {
            p_global[i][j] = 9999.0;
        }
    }

    int **p_total = new int *[num_triangles];
    for (int i = 0; i < num_triangles; i++) {
        p_total[i] = new int [9];
        for (int j = 0; j < 9; j++) {
            p_total[i][j] = 9999;
        }
    }

    int j = 0;
    double x = 0, y = 0;
    for (int i = 0; i < num_triangles; i++) {

        for (int k = 0; k < 3; k++) {
            x = triangles_xy[i][0+2*k];
            y = triangles_xy[i][1+2*k];
            if (!value_comparison(sqrt(x*x+y*y), 1.0, error)) {
                p_total[i][0+3*k] = 0;  // non-boundary
                for (j = 0; j < num1; j++) {
                    if (value_comparison(x, p_internal[j][0], error) && value_comparison(y, p_internal[j][1], error)) {
                        p_total[i][1+3*k] = j;
                        for (int l = 0; l < num3; l++) {
                            if (value_comparison(x, p_global[l][0], error) && value_comparison(y, p_global[l][1], error)) {
                                p_total[i][2+3*k] = l;
                            }
                        }
                        break;
                    }
                }
                if (num1 == 0 || j == num1) {
                    p_total[i][1+3*k] = num1;
                    p_total[i][2+3*k] = num3;
                    p_internal[num1][0] = x;
                    p_internal[num1][1] = y;
                    p_global[num3][0] = x;
                    p_global[num3][1] = y;
                    num1++;
                    num3++;
                }
            }
            else {
                p_total[i][0+3*k] = 1;
                for (j = 0; j < num2; j++) {
                    if (value_comparison(x, p_boundary[j][0], error) && value_comparison(y, p_boundary[j][1], error)) {
                        p_total[i][1+3*k] = j;
                        for (int l = 0; l < num3; l++) {
                            if (value_comparison(x, p_global[l][0], error) && value_comparison(y, p_global[l][1], error)) {
                                p_total[i][2+3*k] = l;
                            }
                        }
                        break;
                    }   
                }
                if (num2 == 0 || j == num2) {
                    p_total[i][1+3*k] = num2;
                    p_total[i][2+3*k] = num3;
                    p_boundary[num2][0] = x;
                    p_boundary[num2][1] = y;
                    p_global[num3][0] = x;
                    p_global[num3][1] = y;
                    num2++;
                    num3++;
                }
            }
        }

    }

    p_total[0] = enlarge_1d(p_total[0], 9, 11);
    p_total[0][9] = num_edges;
    p_total[0][10] = num_vertices;

    for (int i = 0; i < num_vertices-num_edges; i++) {
        delete [] p_internal[i];
    }
    delete [] p_internal;

    for (int i = 0; i < num_edges; i++) {
        delete [] p_boundary[i];
    }
    delete [] p_boundary;

    for (int i = 0; i < num_vertices; i++) {
        delete [] p_global[i];
    }
    delete [] p_global;

    // for (int i = 0; i < num_triangles; i++) {
    //     printf("[%d] %d %d [%d] %d %d [%d] %d %d \n", \
    //     p_total[i][0], p_total[i][1], p_total[i][2], p_total[i][3], p_total[i][4], p_total[i][5], p_total[i][6], p_total[i][7], p_total[i][8]);
    // }
    // // exit(0);

    return p_total;

}

double **triangles_a_b_c_theta_x0_y0(double **triangles_xy) {
    printf("constructing infomation in computational domains...\n");
    int num_triangles = (int) triangles_xy[0][6];

    double *temp1, dx_dy[6];
    double **triangles_uv = new double *[num_triangles];
    for (int j = 0; j < num_triangles; j++) {
        triangles_uv[j] = new double [6];
        for (int k = 0; k < 6; k++) {
            triangles_uv[j][k] = 0;
        }
        triangles_uv[j][3] = atan2(triangles_xy[j][3] - triangles_xy[j][1], triangles_xy[j][2] - triangles_xy[j][0]);
            
        temp1 = foot_of_perpendicular(triangles_xy[j][0], triangles_xy[j][1], triangles_xy[j][2], triangles_xy[j][3], triangles_xy[j][4], triangles_xy[j][5]);
        triangles_uv[j][4] = temp1[0];
        triangles_uv[j][5] = temp1[1];

        dx_dy[0] = temp1[0]-triangles_xy[j][0];
        dx_dy[1] = temp1[1]-triangles_xy[j][1];
        dx_dy[2] = temp1[0]-triangles_xy[j][2];
        dx_dy[3] = temp1[1]-triangles_xy[j][3];
        dx_dy[4] = temp1[0]-triangles_xy[j][4];
        dx_dy[5] = temp1[1]-triangles_xy[j][5];

        triangles_uv[j][0] = sqrt(dx_dy[2] * dx_dy[2] + dx_dy[3] * dx_dy[3]);
        triangles_uv[j][1] = sqrt(dx_dy[0] * dx_dy[0] + dx_dy[1] * dx_dy[1]);
        triangles_uv[j][2] = sqrt(dx_dy[4] * dx_dy[4] + dx_dy[5] * dx_dy[5]);
            
        // delete temp1;
    }

    triangles_uv[0] = enlarge_1d(triangles_uv[0], 6, 7);
    triangles_uv[0][6] = (int) num_triangles;

    // double x1, x2, x3, y1, y2, y3;
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     x1 = triangles_uv[i][4] + (-triangles_uv[i][1]) * cos(triangles_uv[i][3]) - (0) * sin(triangles_uv[i][3]);
    //     y1 = triangles_uv[i][5] + (-triangles_uv[i][1]) * sin(triangles_uv[i][3]) + (0) * cos(triangles_uv[i][3]);

    //     x2 = triangles_uv[i][4] + (triangles_uv[i][0]) * cos(triangles_uv[i][3]) - (0) * sin(triangles_uv[i][3]);
    //     y2 = triangles_uv[i][5] + (triangles_uv[i][0]) * sin(triangles_uv[i][3]) + (0) * cos(triangles_uv[i][3]);

    //     x3 = triangles_uv[i][4] + (0) * cos(triangles_uv[0][3]) - (triangles_uv[i][2]) * sin(triangles_uv[i][3]);
    //     y3 = triangles_uv[i][5] + (0) * sin(triangles_uv[0][3]) + (triangles_uv[i][2]) * cos(triangles_uv[i][3]);

    //     printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n", x1, y1, x2, y2, x3, y3);
    // }
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     printf("a, b, c, theta, x0, y0: %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f;\n", triangles_uv[i][0], triangles_uv[i][1], triangles_uv[i][2], triangles_uv[i][3], triangles_uv[i][4], triangles_uv[i][5]);
    // }
    // exit(0);

    return triangles_uv;
}

double **normal_vector(double **fc_zeta, double **fc_theta, double **triangles_xy, int **points_sequence, int num_sections, int num_modes, int nfp, int nfp_ft) {
    printf("constructing normal vectors on the boundary...\n");

    int num_modesD = 2 * num_modes - 1;
    int num_edges = points_sequence[0][9];
    int num_triangles_in_plane = (int) triangles_xy[0][6];

    double R_zeta, R_theta, Z_zeta, Z_theta;
    double dzeta = 2 * pi / num_sections;
    double *zeta = new double [num_sections]();
    double *theta = new double [num_edges]();;

    double **vector_3_fixed = new double *[3];
    double ***vector_3_sec_edges = new double **[3];
    for (int i = 0; i < 3; i++) {
        vector_3_fixed[i] = new double [num_edges*num_modesD]();
        vector_3_sec_edges[i] = new double *[num_sections];
        for (int j = 0; j < num_sections; j++) {
            vector_3_sec_edges[i][j] = new double [num_edges]();
        }
    }

    for (int i = 0; i < num_sections; i++) {
        zeta[i] = dzeta * i;
    }

    for (int i = 0; i < num_triangles_in_plane; i++) {
        for(int j = 0; j < 3; j++) {
            if (points_sequence[i][0+3*j] == 1) {
                theta[points_sequence[i][1+3*j]] = atan2(triangles_xy[i][1+2*j], triangles_xy[i][0+2*j]);
                if (theta[points_sequence[i][1+3*j]] < 0) {
                    theta[points_sequence[i][1+3*j]] += 2 * pi;
                }
            }
        }
    }

    // for (int i = 0; i < num_edges; i++) {
    //     printf("i = %d, %0.9f\n", i, theta[i]);
    // }
    // exit(0);
    
    double temp1[4];
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                if (points_sequence[j][0+3*k] == 1) {
                    R_zeta =  boundary_R(theta[points_sequence[j][1+3*k]], zeta[i]*nfp, fc_zeta);
                    R_theta = boundary_R(theta[points_sequence[j][1+3*k]], zeta[i]*nfp, fc_theta);
                    Z_zeta =  boundary_Z(theta[points_sequence[j][1+3*k]], zeta[i]*nfp, fc_zeta);
                    Z_theta = boundary_Z(theta[points_sequence[j][1+3*k]], zeta[i]*nfp, fc_theta);

                    temp1[0] = Z_theta;
                    temp1[1] = - R_theta;
                    temp1[2] = R_theta * Z_zeta - R_zeta * Z_theta;

                    temp1[3] = 1.0 / sqrt(temp1[0]*temp1[0] + temp1[1]*temp1[1] + temp1[2]*temp1[2]);
                    vector_3_sec_edges[0][i][points_sequence[j][1+3*k]] = temp1[0] * temp1[3];
                    vector_3_sec_edges[1][i][points_sequence[j][1+3*k]] = temp1[1] * temp1[3];
                    vector_3_sec_edges[2][i][points_sequence[j][1+3*k]] = temp1[2] * temp1[3];

                }
            }
        }
    }

    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < num_edges; j++) {
    //         for (int k = 0 ; k < num_sections; k++) {
    //             printf("%0.9f  ", vector_3_sec_edges[i][k][j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n\n\n");
    // }
    // exit(0);

    for (int i = 0; i < num_edges*num_modesD; i++) {
        for (int j = 0; j < num_sections; j++) {
            vector_3_fixed[0][i] += vector_3_sec_edges[0][j][i/num_modesD] * fourier_basis(i%num_modesD, nfp_ft, zeta[j]) * dzeta;
            vector_3_fixed[1][i] += vector_3_sec_edges[1][j][i/num_modesD] * fourier_basis(i%num_modesD, nfp_ft, zeta[j]) * dzeta;
            vector_3_fixed[2][i] += vector_3_sec_edges[2][j][i/num_modesD] * fourier_basis(i%num_modesD, nfp_ft, zeta[j]) * dzeta;
        }
    }

    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < num_edges*num_modesD; j++) {
    //         printf("%0.9f ", vector_3_fixed[i][j]);
    //         if ((j+1) % num_modesD == 0) {
    //             printf("\n");
    //         }
    //     }
    //     printf("\n\n\n");
    // }
    // exit(0);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < num_sections; j++) {
            delete [] vector_3_sec_edges[i][j];
        }
        delete [] vector_3_sec_edges[i];
    }
    delete [] vector_3_sec_edges;
    delete [] zeta;
    delete [] theta;

    return vector_3_fixed;

}


double *foot_of_perpendicular(double x1, double y1, double x2, double y2, double x3, double y3) {
    double a = y2 - y1, \
           b = x1 - x2, \
           c = x2 * y1 - x1 * y2;
    double * foot = new double [2];
    foot[0] = (b * b * x3 - a * b * y3 - a * c) / (a * a + b * b);
    foot[1] = (a * a * y3 - a * b * x3 - b * c) / (a * a + b * b);
    
    return foot;
}

double ***G_tri_20_18(double **triangles_uv) {
    printf("calculating G...\n");
    
    int num_triangles_in_plane = (int) triangles_uv[0][6];
 
    MatrixXd T1(20, 20);
    // MatrixXd T1_inverse(20, 20);
    MatrixXd T2(20, 18);
    MatrixXd R1(6, 6);
    MatrixXd R(18, 18);
    MatrixXd G(20, 18);
    MatrixXd zeros(6,6);
    MatrixXd E = MatrixXd::Identity(20, 20);
    T1.setZero();
    // T1_inverse.setZero();
    T2.setZero();
    R1.setZero();
    R.setZero();
    G.setZero();
    zeros.setZero();
     
    double a, b, c, cosine, sine;
    double ***g= new double **[num_triangles_in_plane];
    for (int j = 0; j < num_triangles_in_plane; j++) {
        g[j] = new double *[20];
        for (int k = 0; k < 20; k++) {
            g[j][k] = new double [18];
        }
    } 
   
    for (int j = 0; j < num_triangles_in_plane; j++) { 
        a = triangles_uv[j][0];
        b = triangles_uv[j][1];
        c = triangles_uv[j][2];
        cosine = cos(triangles_uv[j][3]);
        sine = sin(triangles_uv[j][3]);
        T1 << 1, -b, 0, b*b, 0, 0, -b*b*b, 0, 0, 0, b*b*b*b, 0, 0, 0, 0, -b*b*b*b*b, 0, 0, 0, 0,
              0, 1, 0, -2*b, 0, 0, 3*b*b, 0, 0, 0, -4*b*b*b, 0, 0, 0, 0, 5*b*b*b*b, 0, 0, 0, 0,
              0, 0, 1, 0, -b, 0, 0, b*b, 0, 0, 0, -b*b*b, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 2, 0, 0, -6*b, 0, 0, 0, 12*b*b, 0, 0, 0, 0, -20*b*b*b, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, -2*b, 0, 0, 0, 3*b*b, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 2, 0, 0, -2*b, 0, 0, 0, 2*b*b, 0, 0, 0, -2*b*b*b, 0, 0, 0,
              1, a, 0, a*a, 0, 0, a*a*a, 0, 0, 0, a*a*a*a, 0, 0, 0, 0, a*a*a*a*a, 0, 0, 0, 0,
              0, 1, 0, 2*a, 0, 0, 3*a*a, 0, 0, 0, 4*a*a*a, 0, 0, 0, 0, 5*a*a*a*a, 0, 0, 0, 0,
              0, 0, 1, 0, a, 0, 0, a*a, 0, 0, 0, a*a*a, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 2, 0, 0, 6*a, 0, 0, 0, 12*a*a, 0, 0, 0, 0, 20*a*a*a, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 2*a, 0, 0, 0, 3*a*a, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 2, 0, 0, 2*a, 0, 0, 0, 2*a*a, 0, 0, 0, 2*a*a*a, 0, 0, 0,
              1, 0, c, 0, 0, c*c, 0, 0, 0, c*c*c, 0, 0, 0, 0, c*c*c*c, 0, 0 ,0 ,0, c*c*c*c*c,
              0, 1, 0, 0, c, 0, 0, 0, c*c, 0, 0, 0, 0, c*c*c, 0, 0, 0, 0, c*c*c*c, 0,
              0, 0, 1, 0, 0, 2*c, 0, 0, 0, 3*c*c, 0, 0, 0, 0, 4*c*c*c, 0, 0, 0, 0, 5*c*c*c*c,
              0, 0, 0, 2, 0, 0, 0, 2*c, 0, 0, 0, 0, 2*c*c, 0, 0, 0, 0, 2*c*c*c, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 0, 2*c, 0, 0, 0, 0, 3*c*c, 0, 0, 0, 0, 4*c*c*c, 0,
              0, 0, 0, 0, 0, 2, 0, 0, 0, 6*c, 0, 0, 0, 0, 12*c*c, 0, 0, 0, 0, 20*c*c*c,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5*a*a*a*a*c, 3*a*a*c*c*c-2*a*a*a*a*c, -2*a*c*c*c*c+3*a*a*a*c*c, c*c*c*c*c-4*a*a*c*c*c, 5*a*c*c*c*c,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5*b*b*b*b*c, 3*b*b*c*c*c-2*b*b*b*b*c, 2*b*c*c*c*c-3*b*b*b*c*c, c*c*c*c*c-4*b*b*c*c*c, -5*b*c*c*c*c;
         
        // T1_inverse = T1.inverse();   

        T2 = (T1.fullPivHouseholderQr().solve(E)).eval().block<20,18>(0,0);      
 
        R1 << 1, 0, 0, 0, 0, 0,
                0, cosine, sine, 0, 0, 0,
                0, -sine, cosine, 0, 0, 0,
                0, 0, 0, cosine*cosine, 2*sine*cosine, sine*sine,
                0, 0, 0, -sine*cosine, cosine*cosine-sine*sine, sine*cosine,
                0, 0, 0, sine*sine, -2*sine*cosine, cosine*cosine;
 
        R << R1, zeros, zeros, zeros, R1, zeros, zeros, zeros, R1;  
         
        G = T2 * R;
        // printf("determinant1: %0.60f\n", T1.determinant());
        // printf("determinant2: %0.60f\n", -64.0*pow(a+b,17)*pow(c,20)*(a*a+c*c)*(b*b+c*c));
        // std::cout << T1 << "\n\n" << T1_inverse << "\n\n" << T1 * T1_inverse << "\n\n" << T2 << "\n\n" << R1 << "\n\n" << R << "\n\n" << G << "\n\n"<< std::endl;
        // exit(0);
         
        for (int k = 0; k < 20; k++) {
            for (int l = 0; l < 18; l++) {
                g[j][k][l] = G(k,l);
            }
        }     
    }

    return g;
}

double ***M_tri_18_18(double ***G, double **triangles_uv) {
    printf("calculating M...\n");
    extern double integral2(double a, double b, double c, double ***G, int i, int j, int num); 
    int num_triangles = (int) triangles_uv[0][6];
    double ***M = new double **[num_triangles];
    for (int i = 0; i < num_triangles; i++) {
        M[i] = new double *[18];
        for (int j = 0; j < 18; j++) {
            M[i][j] = new double [18];
        }
    }
    
    #pragma omp parallel for num_threads (NT)
    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 18; j++) {
            for (int k = 0; k < 18; k++) {
                M[i][j][k] = integral2(triangles_uv[i][0], triangles_uv[i][1], triangles_uv[i][2], G, j, k, i);
            }
        }
    }

    // int i = 0;
    // for (int j = 0; j < 18; j++) {
    //     for (int k = 0; k < 18; k++) {
    //         printf("%0.9f ",M[i][j][k]);
    //     }
    //     printf("\n");
    // }   
    // exit(0);

    return M;
}

double *coordinate_funRZ(double u, double v, double x0, double y0, double zeta, double **fc, double **fc_zeta, double **fc_zetazeta, int num_variables, int nfp, double sine, double cosine) {
    // R  Z  Rzeta Zzeta  Rzetazeta Zzetazeta 
    int num = num_variables;
    double *ptr = new double [num];
    ptr[0] = coordinate_R(u, v, x0, y0, zeta, fc, nfp, sine, cosine);
    ptr[1] = coordinate_Z(u, v, x0, y0, zeta, fc, nfp, sine, cosine);
    // ptr[2] = coordinate_R(u, v, x0, y0, zeta, fc_zeta, nfp, sine, cosine);
    // ptr[3] = coordinate_Z(u, v, x0, y0, zeta, fc_zeta, nfp, sine, cosine);
    // ptr[4] = coordinate_R(u, v, x0, y0, zeta, fc_zetazeta, nfp, sine, cosine);
    // ptr[5] = coordinate_Z(u, v, x0, y0, zeta, fc_zetazeta, nfp, sine, cosine);
    return ptr;
}

double jacobian(double u, double v, double ***R, double ***Z, double ***G, int position_sec, int position_tri) {
    extern int m_order[20];
    extern int n_order[20];

    // double uv[20];
    double uv_u[20], uv_v[20];
    double R_u = 0, R_v = 0;
    double Z_u = 0, Z_v = 0;

    for (int i = 0; i < 20; i++) {
        // uv[i] = power_nonnegative(u, m_order[i]) * power_nonnegative(v, n_order[i]);
        uv_u[i] = m_order[i] * power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] = n_order[i] * power_nonnegative(u, m_order[i]) * power_nonnegative(v, n_order[i]-1);
        for (int j = 0; j < 18; j++) {
            R_u += G[position_tri][i][j] * uv_u[i] * R[position_sec][position_tri][j];
            R_v += G[position_tri][i][j] * uv_v[i] * R[position_sec][position_tri][j];

            Z_u += G[position_tri][i][j] * uv_u[i] * Z[position_sec][position_tri][j];
            Z_v += G[position_tri][i][j] * uv_v[i] * Z[position_sec][position_tri][j];            
        }
    }

    return R_u * Z_v - R_v * Z_u;
}

double **P_9_9(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri) {
    extern int m_order[20];
    extern int n_order[20];

    int dim = 9;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim];
        for (int j = 0; j < dim; j++) {
            p[i][j] = 0;
        }
    }

    double D, E, H, alpha[5], beta[5], gamma[5];
    double A_31, A_32, A_41, A_42, A_51, A_52;
    double uv[20], uv_u[20], uv_v[20], uv_uu[20], uv_uv[20], uv_vv[20];
    double R_u = 0, R_v = 0, R_uu = 0, R_uv = 0, R_vv = 0, R_zetau = 0, R_zetav = 0, R_zetazeta = 0, R_zeta = 0;
    double Z_u = 0, Z_v = 0, Z_uu = 0, Z_uv = 0, Z_vv = 0, Z_zetau = 0, Z_zetav = 0, Z_zetazeta = 0, Z_zeta = 0;

    for (int i = 0; i < 20; i++) {
        uv[i] =                                   power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] *                   power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] *                   power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);
        uv_uu[i] = m_order[i] * (m_order[i]-1) *  power_nonnegative(u, m_order[i]-2) * power_nonnegative(v, n_order[i]);
        uv_uv[i] = m_order[i] *  n_order[i] *     power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]-1);
        uv_vv[i] = n_order[i] * (n_order[i]-1) *  power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-2);
        for (int j = 0; j < 18; j++) {
            R_u +=        G[position_tri][i][j] * uv_u[i] *  R[position_sec][position_tri][j];
            R_v +=        G[position_tri][i][j] * uv_v[i] *  R[position_sec][position_tri][j];
            R_uu +=       G[position_tri][i][j] * uv_uu[i] * R[position_sec][position_tri][j];
            R_uv +=       G[position_tri][i][j] * uv_uv[i] * R[position_sec][position_tri][j];
            R_vv +=       G[position_tri][i][j] * uv_vv[i] * R[position_sec][position_tri][j];

            Z_u +=        G[position_tri][i][j] * uv_u[i] *  Z[position_sec][position_tri][j];
            Z_v +=        G[position_tri][i][j] * uv_v[i] *  Z[position_sec][position_tri][j];
            Z_uu +=       G[position_tri][i][j] * uv_uu[i] * Z[position_sec][position_tri][j];
            Z_uv +=       G[position_tri][i][j] * uv_uv[i] * Z[position_sec][position_tri][j];
            Z_vv +=       G[position_tri][i][j] * uv_vv[i] * Z[position_sec][position_tri][j];

            R_zeta +=     G[position_tri][i][j] * uv[i] *    Rzeta[position_sec][position_tri][j];
            R_zetau +=    G[position_tri][i][j] * uv_u[i] *  Rzeta[position_sec][position_tri][j];
            R_zetav +=    G[position_tri][i][j] * uv_v[i] *  Rzeta[position_sec][position_tri][j];

            Z_zeta +=     G[position_tri][i][j] * uv[i] *    Zzeta[position_sec][position_tri][j];
            Z_zetau +=    G[position_tri][i][j] * uv_u[i] *  Zzeta[position_sec][position_tri][j];
            Z_zetav +=    G[position_tri][i][j] * uv_v[i] *  Zzeta[position_sec][position_tri][j];

            R_zetazeta += G[position_tri][i][j] * uv[i] *    Rzetazeta[position_sec][position_tri][j];

            Z_zetazeta += G[position_tri][i][j] * uv[i] *    Zzetazeta[position_sec][position_tri][j];
            
        }
    }

    D = R_u * Z_v - R_v * Z_u;
    E = R_v * Z_zeta - R_zeta * Z_v;
    H = R_zeta * Z_u - R_u * Z_zeta;

    alpha[0] = (R_zetav * Z_u - R_zetau * Z_v) / D;
    alpha[1] = (Z_zetav * Z_u - Z_zetau * Z_v) / D;
    alpha[2] = - R_zeta;
    alpha[3] = - Z_zeta;
    alpha[4] = 0;

    beta[0] = (R_zetau * R_v - R_zetav * R_u) / D;
    beta[1] = (Z_zetau * R_v - Z_zetav * R_u) / D;
    beta[2] = 0;
    beta[3] = - R_zeta;
    beta[4] = - Z_zeta;

    gamma[0] = - 2 * E/D * R_zetau - 2 * H/D * R_zetav - R_zetazeta;
    gamma[1] = - 2 * E/D * Z_zetau - 2 * H/D * Z_zetav - Z_zetazeta;
    gamma[2] = R_zeta * R_zeta;
    gamma[3] = 2 * R_zeta * Z_zeta;
    gamma[4] = Z_zeta * Z_zeta;

    A_31 = - R_vv*Z_u*Z_u*Z_v + R_v*Z_vv*Z_u*Z_u + 2*R_uv*Z_u*Z_v*Z_v - 2*R_v*Z_uv*Z_u*Z_v - R_uu*Z_v*Z_v*Z_v + R_v*Z_uu*Z_v*Z_v;
    A_32 =   R_vv*Z_u*Z_u*Z_u - 2*R_uv*Z_u*Z_u*Z_v - R_u*Z_vv*Z_u*Z_u + R_uu*Z_u*Z_v*Z_v + 2*R_u*Z_uv*Z_u*Z_v - R_u*Z_uu*Z_v*Z_v;
    A_41 =   R_uu*R_v*Z_v*Z_v - R_u*R_uv*Z_v*Z_v + R_v*R_v*Z_u*Z_uv - R_v*R_v*Z_uu*Z_v - R_u*R_v*Z_u*Z_vv + R_u*R_v*Z_uv*Z_v + R_u*R_vv*Z_u*Z_v - R_uv*R_v*Z_u*Z_v;
    A_42 =   R_uv*R_v*Z_u*Z_u - R_u*R_vv*Z_u*Z_u + R_u*R_u*Z_u*Z_vv - R_u*R_u*Z_uv*Z_v + R_u*R_uv*Z_u*Z_v - R_u*R_v*Z_u*Z_uv + R_u*R_v*Z_uu*Z_v - R_uu*R_v*Z_u*Z_v;
    A_51 =   Z_vv*R_u*R_u*R_v - R_vv*Z_v*R_u*R_u - 2*Z_uv*R_u*R_v*R_v + 2*R_uv*Z_v*R_u*R_v + Z_uu*R_v*R_v*R_v - R_uu*Z_v*R_v*R_v;
    A_52 = - Z_vv*R_u*R_u*R_u + 2*Z_uv*R_u*R_u*R_v + R_vv*Z_u*R_u*R_u - Z_uu*R_u*R_v*R_v - 2*R_uv*Z_u*R_u*R_v + R_uu*Z_u*R_v*R_v;

    p[0][0] = Z_v / D;
    p[0][1] = - Z_u / D;

    p[1][0] = - R_v / D;
    p[1][1] = R_u / D;

    p[2][0] = A_31 / (D*D*D);
    p[2][1] = A_32 / (D*D*D);
    p[2][2] = (Z_v * Z_v) / (D*D);
    p[2][3] = - (2 * Z_u * Z_v) / (D*D);
    p[2][4] = (Z_u * Z_u) / (D*D);

    p[3][0] = A_41 / (D*D*D);
    p[3][1] = A_42 / (D*D*D);
    p[3][2] = - (R_v * Z_v) / (D*D);
    p[3][3] = (R_u * Z_v + R_v * Z_u) / (D*D);
    p[3][4] = - (R_u * Z_u) / (D*D);

    p[4][0] = A_51 / (D*D*D);
    p[4][1] = A_52 / (D*D*D);
    p[4][2] = (R_v * R_v) / (D*D);
    p[4][3] = - (2 * R_u * R_v) / (D*D);
    p[4][4] = (R_u * R_u) / (D*D);

    for (int i = 0; i < 5; i ++) {
        p[5][0] += alpha[i] * p[i][0];
        p[5][1] += alpha[i] * p[i][1];
        p[5][2] += alpha[i] * p[i][2];
        p[5][3] += alpha[i] * p[i][3];
        p[5][4] += alpha[i] * p[i][4];

        p[6][0] += beta[i] * p[i][0];
        p[6][1] += beta[i] * p[i][1];
        p[6][2] += beta[i] * p[i][2];
        p[6][3] += beta[i] * p[i][3];
        p[6][4] += beta[i] * p[i][4];

        p[7][0] += gamma[i] * p[i][0];
        p[7][1] += gamma[i] * p[i][1];
        p[7][2] += gamma[i] * p[i][2];
        p[7][3] += gamma[i] * p[i][3];
        p[7][4] += gamma[i] * p[i][4];
    }

    
    p[5][5] = Z_v / D;
    p[5][6] = - Z_u / D;

   
    p[6][5] = - R_v / D;
    p[6][6] = R_u / D;

    
    p[7][5] = 2 * E / D;
    p[7][6] = 2 * H / D;
    p[7][7] = 1;

    p[8][0] = E / D;
    p[8][1] = H / D;
    p[8][8] = 1;

    return p;
}

double ***F_sec_tri_18(double ***M_ij, double ***f_i, int num_sections, int num_triangles_in_plane) {
    printf("calculating F_1...\n");
 
    MatrixXd x_matrix(18, 1);
    MatrixXd A_matrix(18, 18);
    MatrixXd b_matrix(18, 1);

    double ***F = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        F[i] = new double *[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            F[i][j] = new double [18];

            for (int k = 0; k < 18; k++) {
                b_matrix(k,0) = f_i[i][j][k];
                for (int l = 0; l < 18; l++) {
                    A_matrix(k,l) = M_ij[j][k][l];
                }
            }

            x_matrix = A_matrix.fullPivHouseholderQr().solve(b_matrix);
            // x_matrix = A_matrix.fullPivLu().solve(b_matrix);
            // x_matrix = A_matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_matrix);

            for (int k = 0; k < 18; k++) {
                F[i][j][k] = x_matrix(k,0);
            }

        }
    } 

    // int num = 0;
    // for (int i = 0; i < num_triangles_in_plane; i++) {
    //     for (int j = 0; j < 18; j++) {
    //         printf("%f ",F[num][i][j]);
    //     }
    //     printf("\n");
    // }

    return F;
}

double ****F_sec_tri_18_num(double ***M_ij, double ****f_i, int num_sections, int num_triangles_in_plane, int num_variables) {
    printf("calculating F_1d...\n");

    int num = num_variables;

    MatrixXd x_matrix(18, 1);
    MatrixXd A_matrix(18, 18);
    MatrixXd b_matrix(18, 1);

    double ****F = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        F[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {

            F[i][j] = new double *[18];
            for (int k = 0; k < 18; k++) {
                F[i][j][k] = new double [num];
                for (int l = 0; l < 18; l++) {
                    A_matrix(k,l) = M_ij[j][k][l];
                }
                for (int l = 0; l < num; l++) {
                    F[i][j][k][l] = 0;
                }
            }

            for (int l = 0; l < num; l++) {
                for (int k = 0; k < 18; k++) {
                    b_matrix(k,0) = f_i[i][j][k][l];
                }
                x_matrix = A_matrix.fullPivHouseholderQr().solve(b_matrix);
                for (int k = 0; k < 18; k++) {
                    F[i][j][k][l] = x_matrix(k,0);
                } 
            }

        }
    }

    return F;
}

double *****F_sec_tri_18_m_m(double ***M_ij, double *****f_i, int num_sections, int num_triangles_in_plane, int rank) {
    printf("calculating F_2d...\n");
 
    MatrixXd x_matrix(18, 1);
    MatrixXd A_matrix(18, 18);
    MatrixXd b_matrix(18, 1);

    double *****F = new double ****[num_sections];
    for (int i = 0; i < num_sections; i++) {
        F[i] = new double ***[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {

            F[i][j] = new double **[18];
            for (int k = 0; k < 18; k++) {
                F[i][j][k] = new double *[rank];
                for (int l = 0; l < 18; l++) {
                    A_matrix(k,l) = M_ij[j][k][l];
                }
                for (int s = 0; s < rank; s++) {
                    F[i][j][k][s] = new double [rank];
                    for (int t = 0; t < rank; t++) {
                        F[i][j][k][s][t] = 0;   
                    }
                }    
            }

            for (int s = 0; s < rank; s++) {
                for (int t = 0; t < rank; t++) {
                    for (int k = 0; k < 18; k++) {
                        b_matrix(k,0) = f_i[i][j][k][s][t];
                    }
                    x_matrix = A_matrix.fullPivHouseholderQr().solve(b_matrix);
                    for (int k = 0; k < 18; k++) {
                        F[i][j][k][s][t] = x_matrix(k,0);
                    } 
                }
            }

        }
    }

    return F;
}

double ******F_sec_tri_18_derivates_m_m(double ***M_ij, double ******f_i, int num_sections, int num_triangles_in_plane, int rank, int num_derivates) {
    printf("calculating F_3d...\n");
 
    MatrixXd x_matrix(18, 1);
    MatrixXd A_matrix(18, 18);
    MatrixXd b_matrix(18, 1);

    double ******F = new double *****[num_sections];
    for (int i = 0; i < num_sections; i++) {
        F[i] = new double ****[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {

            F[i][j] = new double ***[18];
            for (int k = 0; k < 18; k++) {
                F[i][j][k] = new double **[num_derivates];

                for (int l = 0; l < 18; l++) {
                    A_matrix(k,l) = M_ij[j][k][l];
                }

                for (int l = 0; l < num_derivates; l++) {
                    F[i][j][k][l] = new double *[rank];
                    for (int s = 0; s < rank; s++) {
                        F[i][j][k][l][s] = new double [rank]();
                    } 
                }
                   
            }

            for (int s = 0; s < rank; s++) {
                for (int t = 0; t < rank; t++) {
                    for (int l = 0; l < num_derivates; l++) {
                        for (int k = 0; k < 18; k++) {
                            b_matrix(k,0) = f_i[i][j][k][l][s][t];
                        }
                        x_matrix = A_matrix.fullPivHouseholderQr().solve(b_matrix);
                        for (int k = 0; k < 18; k++) {
                            F[i][j][k][l][s][t] = x_matrix(k,0);
                        } 
                    }
                    
                }
            }

        }
    }

    return F;
}









// double **re_axis(double **triangles_xy, double **fourier_coefficients, int nfp, int n, int order) {
//     int num_coefficients = (int) fourier_coefficients[0][6];

//     double zeta[n], theta[n];
//     for (int i = 0; i < n; i++) {
//         zeta[i] = i * (2*pi/n);
//         theta[i] = i * (2*pi/n);
//     }

//     double **boundary_points_r = new double *[n];
//     double **boundary_points_z = new double *[n];
//     double sum_r, sum_z;
//     for (int i = 0; i < n; i++) {
//         //zeta
//         boundary_points_r[i] = new double [n];
//         boundary_points_z[i] = new double [n];      
//         for (int j = 0; j < n; j++) {
//             //theta
//             sum_r = 0;
//             sum_z = 0;
//             for (int k = 0; k < num_coefficients; k++) {
//                 sum_r += fourier_coefficients[k][2] * cos(fourier_coefficients[k][1]*theta[j] - nfp*fourier_coefficients[k][0]*zeta[i]) + \
//                          fourier_coefficients[k][4] * sin(fourier_coefficients[k][1]*theta[j] - nfp*fourier_coefficients[k][0]*zeta[i]);
//                 sum_z += fourier_coefficients[k][3] * sin(fourier_coefficients[k][1]*theta[j] - nfp*fourier_coefficients[k][0]*zeta[i]) + \
//                          fourier_coefficients[k][5] * cos(fourier_coefficients[k][1]*theta[j] - nfp*fourier_coefficients[k][0]*zeta[i]);
//             }
//             boundary_points_r[i][j] = sum_r;
//             boundary_points_z[i][j] = sum_z;
//         }
//     }

//     double *origin_r = new double [n];
//     double *origin_z = new double [n];
//     for (int i = 0; i < n; i++) {
//         sum_r = 0;
//         sum_z = 0;
//         for (int j = 0; j < num_coefficients; j++) {
//             sum_r += 0.5*(fourier_coefficients[j][2] * cos(fourier_coefficients[j][1]*pi - nfp*fourier_coefficients[j][0]*zeta[i]) + \
//                           fourier_coefficients[j][4] * sin(fourier_coefficients[j][1]*pi - nfp*fourier_coefficients[j][0]*zeta[i])) + \
//                      0.5*(fourier_coefficients[j][2] * cos(fourier_coefficients[j][1]*0 - nfp*fourier_coefficients[j][0]*zeta[i]) + \
//                           fourier_coefficients[j][4] * sin(fourier_coefficients[j][1]*0 - nfp*fourier_coefficients[j][0]*zeta[i]));
//             sum_z += 0.5*(fourier_coefficients[j][3] * sin(fourier_coefficients[j][1]*1.5*pi - nfp*fourier_coefficients[j][0]*zeta[i]) + \
//                           fourier_coefficients[j][5] * cos(fourier_coefficients[j][1]*1.5*pi - nfp*fourier_coefficients[j][0]*zeta[i])) + \
//                      0.5*(fourier_coefficients[j][3] * sin(fourier_coefficients[j][1]*0.5*pi - nfp*fourier_coefficients[j][0]*zeta[i]) + \
//                           fourier_coefficients[j][5] * cos(fourier_coefficients[j][1]*0.5*pi - nfp*fourier_coefficients[j][0]*zeta[i]));
//         }
//         origin_r[i] = sum_r;
//         origin_z[i] = sum_z;
//     }

//     double **coefficients = new double *[order+1];
//     double sum1, sum2, sum3, sum4;
//     for (int i = 0; i < order+1; i++) {
//         coefficients[i] = new double [4];
//         sum1 = 0;
//         sum2 = 0;
//         sum3 = 0;
//         sum4 = 0;
//         for (int j = 0; j < n; j++) {
//             sum1 += origin_r[j] * (2*pi/n) * cos(i*zeta[j]) / pi;
//             sum2 += origin_z[j] * (2*pi/n) * sin(i*zeta[j]) / pi;
//             sum3 += origin_r[j] * (2*pi/n) * sin(i*zeta[j]) / pi;
//             sum4 += origin_z[j] * (2*pi/n) * cos(i*zeta[j]) / pi;
//         }
//         coefficients[i][0] = sum1;
//         coefficients[i][1] = sum2;
//         coefficients[i][2] = sum3;
//         coefficients[i][3] = sum4;
//     }
//     coefficients[0][0] /= 2;
//     coefficients[0][1] /= 2;
//     coefficients[0][2] /= 2;
//     coefficients[0][3] /= 2;
//     coefficients[0] = enlarge_1d(coefficients[0], 4, 5);
//     coefficients[0][4] = order;

//     // for (int i = 0; i < order+1; i++) {
//     //     printf("%d %0.9f %0.9f %0.9f %0.9f\n", i, coefficients[i][0], coefficients[i][1], coefficients[i][2] = sum3, coefficients[i][3]);
//     // }

//     // for (int i = 0; i < n; i++) {
//     //     printf("%0.9f %0.9f\n", origin_r[i], origin_z[i]);
//     // }


//     return coefficients;
// }

// double ***triangles_transform2(double **triangles_xy, double **fourier_coefficients, double **axis_coefficients, int num_sections, int nfp) {
//         int num_triangles = (int) triangles_xy[0][6];
//         int num_coefficients = (int) fourier_coefficients[0][6];
//         int num_order = (int) axis_coefficients[0][4];

//         double **angles = new double *[num_triangles];
//         double **rhos = new double *[num_triangles];
//         for (int i = 0; i < num_triangles; i++) {
//             angles[i] = new double [3];
//             rhos[i] = new double [3];
//             for (int j = 0; j < 3; j++) {
//                 angles[i][j] = atan2(triangles_xy[i][2*j+1], triangles_xy[i][2*j+0]);
//                 rhos[i][j] = sqrt(triangles_xy[i][2*j+1] * triangles_xy[i][2*j+1] + triangles_xy[i][2*j+0] * triangles_xy[i][2*j+0]);
//             }
//         }

//         double ***triangles_RZ = new double **[num_sections];
//         double sum1_r, sum1_z, sum2_r, sum2_z;
//         for (int i = 0; i < num_sections; i++) {
//             triangles_RZ[i] = new double *[num_triangles];
//             for (int j = 0; j < num_triangles; j++ ) {
//                 triangles_RZ[i][j] = new double [6];
//                 for (int k = 0; k < 3; k++) {
//                     sum1_r = 0;
//                     sum1_z = 0;
//                     sum2_r = 0;
//                     sum2_z = 0;
//                     for (int l = 0; l < num_coefficients; l++) {
//                         sum1_r += fourier_coefficients[l][2] * cos(fourier_coefficients[l][1]*angles[j][k] - nfp*fourier_coefficients[l][0]*(2*pi)/num_sections*i) + \
//                                  fourier_coefficients[l][4] * sin(fourier_coefficients[l][1]*angles[j][k] - nfp*fourier_coefficients[l][0]*(2*pi)/num_sections*i);
//                         sum1_z += fourier_coefficients[l][3] * sin(fourier_coefficients[l][1]*angles[j][k] - nfp*fourier_coefficients[l][0]*(2*pi)/num_sections*i) + \
//                                  fourier_coefficients[l][5] * cos(fourier_coefficients[l][1]*angles[j][k] - nfp*fourier_coefficients[l][0]*(2*pi)/num_sections*i);
//                     }
//                     for (int l = 0; l < num_order; l++) {
//                         sum2_r += axis_coefficients[l][0] * cos(l*(2*pi)/num_sections*i) + \
//                                   axis_coefficients[l][2] * sin(l*(2*pi)/num_sections*i);
//                         sum2_z += axis_coefficients[l][1] * sin(l*(2*pi)/num_sections*i) + \
//                                   axis_coefficients[l][3] * cos(l*(2*pi)/num_sections*i);
//                     }
                    
//                     triangles_RZ[i][j][2*k+0] = rhos[j][k] * (sum1_r - sum2_r) + sum2_r;
//                     triangles_RZ[i][j][2*k+1] = rhos[j][k] * (sum1_z - sum2_z) + sum2_z;
//                 }
//             }
//         }

//     for (int i = 0; i < num_sections; i++) {
//         delete[] angles[i];
//         delete[] rhos[i];
//     }
//     delete[] angles;
//     delete[] rhos;

//     for (int i = 0; i < num_triangles; i++) {
//         printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n",triangles_RZ[0][i][0],triangles_RZ[0][i][1],triangles_RZ[0][i][2],triangles_RZ[0][i][3],triangles_RZ[0][i][4],triangles_RZ[0][i][5]);
//     }

//     return triangles_RZ ;
// }

// double **boundary_recalculation(double **fourier_coefficients, int nfp, int n) {
//     int num_coefficients = (int) fourier_coefficients[0][6];

//     double zeta[n], theta[n];
//     for (int i = 0; i < n; i++) {
//         zeta[i] = i * (2*pi/n);
//         theta[i] = i * (2*pi/n);
//     }

//     double **boundary_points_r = new double *[n];
//     double **boundary_points_z = new double *[n];
//     double sum_r, sum_z;
//     for (int i = 0; i < n; i++) {
//         //zeta
//         boundary_points_r[i] = new double [n];
//         boundary_points_z[i] = new double [n];      
//         for (int j = 0; j < n; j++) {
//             //theta
//             sum_r = 0;
//             sum_z = 0;
//             for (int k = 0; k < num_coefficients; k++) {
//                 sum_r += fourier_coefficients[k][2] * cos(fourier_coefficients[k][1]*theta[j] - nfp*fourier_coefficients[k][0]*zeta[i]) + \
//                          fourier_coefficients[k][4] * sin(fourier_coefficients[k][1]*theta[j] - nfp*fourier_coefficients[k][0]*zeta[i]);
//                 sum_z += fourier_coefficients[k][3] * sin(fourier_coefficients[k][1]*theta[j] - nfp*fourier_coefficients[k][0]*zeta[i]) + \
//                          fourier_coefficients[k][5] * cos(fourier_coefficients[k][1]*theta[j] - nfp*fourier_coefficients[k][0]*zeta[i]);
//             }
//             boundary_points_r[i][j] = sum_r;
//             boundary_points_z[i][j] = sum_z;
//         }
//     }

//     double *origin_r = new double [n];
//     double *origin_z = new double [n];
//     for (int i = 0; i < n; i++) {
//         sum_r = 0;
//         sum_z = 0;
//         for (int j = 0; j < num_coefficients; j++) {
//             sum_r += 0.5*(fourier_coefficients[j][2] * cos(fourier_coefficients[j][1]*pi - nfp*fourier_coefficients[j][0]*zeta[i]) + \
//                           fourier_coefficients[j][4] * sin(fourier_coefficients[j][1]*pi - nfp*fourier_coefficients[j][0]*zeta[i])) + \
//                      0.5*(fourier_coefficients[j][2] * cos(fourier_coefficients[j][1]*0 - nfp*fourier_coefficients[j][0]*zeta[i]) + \
//                           fourier_coefficients[j][4] * sin(fourier_coefficients[j][1]*0 - nfp*fourier_coefficients[j][0]*zeta[i]));
//             sum_z += 0.5*(fourier_coefficients[j][3] * sin(fourier_coefficients[j][1]*1.5*pi - nfp*fourier_coefficients[j][0]*zeta[i]) + \
//                           fourier_coefficients[j][5] * cos(fourier_coefficients[j][1]*1.5*pi - nfp*fourier_coefficients[j][0]*zeta[i])) + \
//                      0.5*(fourier_coefficients[j][3] * sin(fourier_coefficients[j][1]*0.5*pi - nfp*fourier_coefficients[j][0]*zeta[i]) + \
//                           fourier_coefficients[j][5] * cos(fourier_coefficients[j][1]*0.5*pi - nfp*fourier_coefficients[j][0]*zeta[i]));
//         }
//         origin_r[i] = sum_r;
//         origin_z[i] = sum_z;
//     }
 
//     double **theta_re = new double *[n];
//     for (int i = 0; i < n; i++) {
//         theta_re[i] = new double [n];        
//         for (int j = 0; j < n; j++) {           
//             theta_re[i][j] = atan2(boundary_points_z[i][j] - origin_z[i], boundary_points_r[i][j] - origin_r[i]);   
//             // theta_re[i][j] = theta[j];       
//             if (theta_re[i][j] < 0) {
//                 theta_re[i][j] += 2*pi;
//             }
//         }
//     }

//     double **fourier_coefficients_re = new double *[num_coefficients];      
//     fourier_coefficients_re[0] = new double [7];    
//     fourier_coefficients_re[0][6] = num_coefficients;  
//     for (int i = 1; i < num_coefficients; i++) {
//         fourier_coefficients_re[i] = new double [6];
//     }
//     double sum1[n], sum2[n], sum3[n], sum4[n];
//     double ssum1, ssum2, ssum3, ssum4;
//     double dtheta, dzeta, cosineang, sineang;
//     for (int i = 0; i < num_coefficients; i++) {
//         for (int j = 0; j < 6; j++) {
//             fourier_coefficients_re[i][j] = fourier_coefficients[i][j];
//         }    

//         ssum1 = 0;
//         ssum2 = 0;
//         ssum3 = 0;
//         ssum4 = 0;

//         for (int k = 0; k < n; k++) {

//             sum1[k] = 0;
//             sum2[k] = 0;
//             sum3[k] = 0;
//             sum4[k] = 0;

//             if (k != n - 1 ) {
//                 dzeta = zeta[k+1] - zeta[k];
//             }
//             else {
//                 dzeta = 2*pi - zeta[k];
//             }

//             for (int l = 0; l < n; l++) {

//                 cosineang = cos(fourier_coefficients_re[i][1] * theta_re[k][l] - nfp*fourier_coefficients_re[i][0] * zeta[k]);
//                 sineang = sin(fourier_coefficients_re[i][1] * theta_re[k][l] - nfp*fourier_coefficients_re[i][0] * zeta[k]);

//                 if (l != n - 1 ) {
//                     dtheta = theta_re[k][l+1] - theta_re[k][l];
//                 }
//                 else {
//                     dtheta = 2*pi - theta_re[k][l];
//                 } 

//                 if (fourier_coefficients_re[i][1] == 0 && fourier_coefficients_re[i][0] == 0) {
//                     sum1[k] += boundary_points_r[k][l] * dtheta * cosineang / (4*pi*pi);
//                     sum2[k] += boundary_points_r[k][l] * dtheta * sineang / (4*pi*pi);
//                     sum3[k] += boundary_points_z[k][l] * dtheta * sineang / (4*pi*pi);
//                     sum4[k] += boundary_points_z[k][l] * dtheta * cosineang / (4*pi*pi);
//                 }
//                 else  {
//                     sum1[k] += boundary_points_r[k][l] * dtheta * cosineang / (2*pi*pi);
//                     sum2[k] += boundary_points_r[k][l] * dtheta * sineang / (2*pi*pi);
//                     sum3[k] += boundary_points_z[k][l] * dtheta * sineang / (2*pi*pi);
//                     sum4[k] += boundary_points_z[k][l] * dtheta * cosineang / (2*pi*pi);
//                 }
                    
//             }

//             ssum1 += sum1[k] * dzeta;
//             ssum2 += sum2[k] * dzeta;
//             ssum3 += sum3[k] * dzeta;
//             ssum4 += sum4[k] * dzeta;
                
//         }

//         fourier_coefficients_re[i][2] = ssum1;
//         fourier_coefficients_re[i][3] = ssum3;
//         fourier_coefficients_re[i][4] = ssum2;
//         fourier_coefficients_re[i][5] = ssum4;
        
//     }
    

//     for (int i = 0; i < (int) fourier_coefficients_re[0][6]; i++) {
//         // printf("%d %d %0.9f %0.9f %0.9f %0.9f\n",(int) fourier_coefficients_re[i][0], (int) fourier_coefficients_re[i][1],\
//         // fourier_coefficients_re[i][2]-fourier_coefficients[i][2], \
//         // fourier_coefficients_re[i][3]-fourier_coefficients[i][3], \
//         // fourier_coefficients_re[i][4]-fourier_coefficients[i][4], \
//         // fourier_coefficients_re[i][5]-fourier_coefficients[i][5]);
//         printf("%d %d %0.9f %0.9f %0.9f %0.9f\n",(int) fourier_coefficients_re[i][0], (int) fourier_coefficients_re[i][1],\
//         fourier_coefficients_re[i][2], \
//         fourier_coefficients_re[i][3], \
//         fourier_coefficients_re[i][4], \
//         fourier_coefficients_re[i][5]);
//     }


//     // for (int i = 0; i < n; i++) {
//     //     printf("%0.9f %0.9f \n",boundary_points_r[0][i],boundary_points_z[0][i]);
//     // }
//     // for (int i = 0; i < n; i++) {
//     //     printf("%0.9f %0.9f \n",boundary_points_r[45][i],boundary_points_z[45][i]);
//     // }
//     // for (int i = 0; i < n; i++) {
//     //     printf("%0.9f %0.9f \n",boundary_points_r[90][i],boundary_points_z[90][i]);
//     // }
//     // for (int i = 0; i < n; i++) {
//     //     printf("%0.9f %0.9f \n",boundary_points_r[135][i],boundary_points_z[135][i]);
//     // }        
//     for (int i = 0; i < n; i++) {
//         delete[] boundary_points_r[i];
//         delete[] boundary_points_z[i];
//         delete[] theta_re[i];
//     }    
//     delete[] boundary_points_r;
//     delete[] boundary_points_z;
//     delete[] theta_re;
//     delete[] origin_r;
//     delete[] origin_z;


//     return fourier_coefficients_re;
// }




// double test1(double u, double v, double theta, double x0, double y0, double zeta, double **fourier_coefficients, int nfp, double sine, double cosine) {
//     return 17.0;
// }
