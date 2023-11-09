#include "input_output.cpp"
#include "triangulation.cpp"
#include "refinement.cpp"
#include "meshes.cpp"
#include "integrals.cpp"
#include "integrals2.cpp"
#include "infrastructure.cpp"
#include "superstructure.cpp"
#include "solvers.cpp"
#include "header.hpp"

int main () {

    const char *file_boundary = "input_boundary.txt";
    const char *file_coils = "coil_cfqs.txt";
    int num_sections = 180;
    int num_segments = 30;
    int nfp = 2;
    double error = 1.0E-14;
    int num_modes = 13;
    int nfp_ft = 2;
    int coil_interpolation = 29;
    int GL_order = 1; // GL_order += 1 for divergence cleaning 
    double beltrami_mu = 0.0; // 0.1 - iota = 1/3
    vector<double> initial_pressure = {10000.0, -10000.0};  // 10000.0, -10000.0

    // // test_FT
    // double temp_ft[50] = {2,2.20938774960796,2.32664351527421,2.34329221923835,2.25646140052488,2.06909050504505,1.78978587537900,1.43232992589165,1.01487455587994,0.558868544592172,0.0877852522924734,-0.374270618617849,-0.803990948220860,-1.18023845440810,-1.48528235315865,-1.70581924104237,-1.83371652408955,-1.86643313901009,-1.80709255506310,-1.66420483632037,-1.45105651629515,-1.18480781056688,-0.885354964395529,-0.574030155186484,-0.272221427350702,-2.44929359829471e-16,0.225158346979008,0.389477193016948,0.483739247461849,0.503848040437148,0.451056516295154,0.331848620536174,0.157481946394274,-0.0567790340780504,-0.292690038537974,-0.530248736457422,-0.749033247789298,-0.929571987279487,-1.05465741534947,-1.11051972398720,-1.08778525229247,-0.982157940959407,-0.794779549052095,-0.532244575565725,-0.206267581477548,0.166977472454741,0.567805549520842,0.974198007380975,1.36313616707078,1.71200797527825};
    // int nummodes = 12;
    // double dZeta = 2 * pi / 50, Zeta = 0;

    // double **test_sec_18 = new double *[50];
    // for (int j = 0; j < 50; j++) {
    //     test_sec_18[j] = new double [18];
    //     for (int k = 0; k < 18; k++) {
    //         test_sec_18[j][k] = temp_ft[j];
    //     }
    // }

    // double *test_18N = new double [18*nummodes];
    // for (int j = 0; j < 18*nummodes; j++) {
    //     test_18N[j] = 0;
    //     for (int l = 0; l < 50; l++) {
    //         Zeta = dZeta * l;
    //         test_18N[j] += test_sec_18[l][j/nummodes] * fourier_basis(j%nummodes, (int) 1, Zeta) *dZeta;
 
    //     }
    // }

    // for (int j = 0; j < 18*nummodes; j++) {
    //     printf("%0.9f  ", test_18N[j]);
    //     if (j%nummodes == nummodes-1) {
    //         printf("\n");
    //     }
    // }
    // printf("\n");
    // exit(0);
    // //


    // //test1
    // int interval = 1000000;
    // double step = 2*pi / interval;
    // double int1 = 0;
    // int Nfp = 2, Ntor = 20;
    // int t = 0*Ntor+4, k = 17*Ntor+17, j = 9*Ntor+13;
    // int Gj = 1;
    // for (int i = 0; i < interval; i++) {
    //     int1 += (sin(Nfp*(t%Ntor)/2.0 * i*step)) * (cos(Nfp*(k%Ntor+1)/2.0 * i*step)) * (- power_nonnegative(Nfp*(j%Ntor+1)/2.0,1) * sin(Nfp*(j%Ntor+1)/2.0 * i*step));

    // }
    // int1 *= step;
    // double int2 = integral3(Nfp, Ntor, t, k, j, Gj, 10e-14);
    // printf("int1 = %0.16f\nint2 = %0.16f\n", int1, int2);
    // exit(0);
    // //
    time_t time_s = time(NULL);
    table();
    mn_order();
    omp_set_num_threads(NT);
    Eigen::initParallel();
    // Eigen::setNbThreads(NT);
    printf("Eigen::nbThreads = %d\n",Eigen::nbThreads());

    // printf("%16.0f\n", factorial(11));
    // printf("%16.0f\n", factorial(12));
    // printf("%16.0f\n", factorial(13));
    // printf("%16.0f\n", factorial(14));
    // printf("%16.0f\n", factorial(15));
    // printf("%16.0f\n", factorial(16));
    // printf("%16.0f\n", factorial(17));
    // printf("%f\n", power_nonnegative(-5.0, 4));
    // printf("%f\n", power_nonnegative(-5.0, 5));
    // printf("%f\n", power_nonnegative(13.1, 0));
    // printf("%f\n", power_nonnegative(-9.2, 0));
    // printf("%f\n", power_nonnegative(1.1, 20));
    // printf("%f\n", power_nonnegative(-1.1, 19));
    // printf("%0.32f %0.32f\n", sin(8765.3), sine_lookup_approximation(8765.3));
    // printf("%0.32f %0.32f\n", cos(8765.3), cosine_lookup_approximation(8765.3));
    // exit(0);

    // double temp_mn[18] = {0.318659991, -0.084175072, -0.182309566, 0.216833459, -0.324608531, 0.201245320, 0.293757175, -0.118708199, -0.162763341, 0.281758244, -0.328692697, 0.251798892, 0.305397449, -0.134624541, -0.140115703, 0.267387792, -0.364668724, 0.255344229};
    // double temp_value = 0;
    // double u = -0.095093462, v = 0.001615232;
    // for (int i = 0; i < 18; i++) {
    //     temp_value += temp_mn[i] * power_nonnegative(u, m_order[i]) * power_nonnegative(v, n_order[i]);
    // }
    // printf("temp_value = %0.9f\n", temp_value);
    // exit(0);

    timing(time_s);
    // extern double sine_table[10000];
    // extern double cosine_table[10000];
    

    // MatrixXd T1(20, 20);
    // double a = 0.0436194, b = 0.0536194, c = 1.83488;
    // T1 << 1, -b, 0, b*b, 0, 0, -b*b*b, 0, 0, 0, b*b*b*b, 0, 0, 0, 0, -b*b*b*b*b, 0, 0, 0, 0,
    //           0, 1, 0, -2*b, 0, 0, 3*b*b, 0, 0, 0, -4*b*b*b, 0, 0, 0, 0, 5*b*b*b*b, 0, 0, 0, 0,
    //           0, 0, 1, 0, -b, 0, 0, b*b, 0, 0, 0, -b*b*b, 0, 0, 0, 0, 0, 0, 0, 0,
    //           0, 0, 0, 2, 0, 0, -6*b, 0, 0, 0, 12*b*b, 0, 0, 0, 0, -20*b*b*b, 0, 0, 0, 0,
    //           0, 0, 0, 0, 1, 0, 0, -2*b, 0, 0, 0, 3*b*b, 0, 0, 0, 0, 0, 0, 0, 0,
    //           0, 0, 0, 0, 0, 2, 0, 0, -2*b, 0, 0, 0, 2*b*b, 0, 0, 0, -2*b*b*b, 0, 0, 0,
    //           1, a, 0, a*a, 0, 0, a*a*a, 0, 0, 0, a*a*a*a, 0, 0, 0, 0, a*a*a*a*a, 0, 0, 0, 0,
    //           0, 1, 0, 2*a, 0, 0, 3*a*a, 0, 0, 0, 4*a*a*a, 0, 0, 0, 0, 5*a*a*a*a, 0, 0, 0, 0,
    //           0, 0, 1, 0, a, 0, 0, a*a, 0, 0, 0, a*a*a, 0, 0, 0, 0, 0, 0, 0, 0,
    //           0, 0, 0, 2, 0, 0, 6*a, 0, 0, 0, 12*a*a, 0, 0, 0, 0, 20*a*a*a, 0, 0, 0, 0,
    //           0, 0, 0, 0, 1, 0, 0, 2*a, 0, 0, 0, 3*a*a, 0, 0, 0, 0, 0, 0, 0, 0,
    //           0, 0, 0, 0, 0, 2, 0, 0, 2*a, 0, 0, 0, 2*a*a, 0, 0, 0, 2*a*a*a, 0, 0, 0,
    //           1, 0, c, 0, 0, c*c, 0, 0, 0, c*c*c, 0, 0, 0, 0, c*c*c*c, 0, 0 ,0 ,0, c*c*c*c*c,
    //           0, 1, 0, 0, c, 0, 0, 0, c*c, 0, 0, 0, 0, c*c*c, 0, 0, 0, 0, c*c*c*c, 0,
    //           0, 0, 1, 0, 0, 2*c, 0, 0, 0, 3*c*c, 0, 0, 0, 0, 4*c*c*c, 0, 0, 0, 0, 5*c*c*c*c,
    //           0, 0, 0, 2, 0, 0, 0, 2*c, 0, 0, 0, 0, 2*c*c, 0, 0, 0, 0, 2*c*c*c, 0, 0,
    //           0, 0, 0, 0, 1, 0, 0, 0, 2*c, 0, 0, 0, 0, 3*c*c, 0, 0, 0, 0, 4*c*c*c, 0,
    //           0, 0, 0, 0, 0, 2, 0, 0, 0, 6*c, 0, 0, 0, 0, 12*c*c, 0, 0, 0, 0, 20*c*c*c,
    //           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5*a*a*a*a*c, 3*a*a*c*c*c-2*a*a*a*a*c, -2*a*c*c*c*c+3*a*a*a*c*c, c*c*c*c*c-4*a*a*c*c*c, 5*a*c*c*c*c,
    //           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5*b*b*b*b*c, 3*b*b*c*c*c-2*b*b*b*b*c, 2*b*c*c*c*c-3*b*b*b*c*c, c*c*c*c*c-4*b*b*c*c*c, -5*b*c*c*c*c;
    // printf("determinant1: %0.32f\n", T1.determinant());
    // printf("determinant2: %0.32f\n", -64*pow(a+b,17)*pow(c,20)*(a*a+c*c)*(b*b+c*c));
    // exit(0);

    if (num_modes % 2 == 0) {
        printf("mode number must be odd\n");
        exit(0);
    }
    
    double **fourier_coefficients = read_boundary(file_boundary);
    double **fc_zeta = fourier_coefficients_dzeta(fourier_coefficients, (int) nfp);
    double **fc_theta = fourier_coefficients_dtheta(fourier_coefficients, (int) nfp);
    double **fc_zetazeta= fourier_coefficients_dzeta(fc_zeta, (int) nfp);
    timing(time_s);

    double ***coils_p = read_coils_p(file_coils, (int) coil_interpolation, 0.5);
    timing(time_s);

    // double *temp1 = biot_savart(1.32502, 0.0, 0.0, coils_p);
    // printf("%0.9f  %0.9f  %0.9f\n", temp1[0], temp1[1], temp1[2]);
    // exit(0);
    
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

    
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     printf("a, b, c = %f, %f, %f\n", triangles_uv[i][0], triangles_uv[i][1], triangles_uv[i][2]);
    // }
    // exit(0);

    // double zeta;
    // double r1, r2, r3, z1, z2, z3;
    // double temp1;
    // double sine, cosine;
    // // for (int j = 0; j < (num_sections*1); j++) {
    // //     zeta = 2*pi / (num_sections*1) * j;
    //     zeta = pi;
    //     for (int i = 0; i < (int) triangles_uv[0][6]; i++) {

    //         sine = sin(triangles_uv[i][3]);
    //         cosine = cos(triangles_uv[i][3]);

    //         r1 = coordinate_R(-triangles_uv[i][1], 0, triangles_uv[i][3], triangles_uv[i][4], triangles_uv[i][5], zeta, fourier_coefficients, nfp, sine, cosine);
    //         r2 = coordinate_R( triangles_uv[i][0], 0, triangles_uv[i][3], triangles_uv[i][4], triangles_uv[i][5], zeta, fourier_coefficients, nfp, sine, cosine);
    //         r3 = coordinate_R( 0, triangles_uv[i][2], triangles_uv[i][3], triangles_uv[i][4], triangles_uv[i][5], zeta, fourier_coefficients, nfp, sine, cosine);

    //         z1 = coordinate_Z(-triangles_uv[i][1], 0, triangles_uv[i][3], triangles_uv[i][4], triangles_uv[i][5], zeta, fourier_coefficients, nfp, sine, cosine);
    //         z2 = coordinate_Z( triangles_uv[i][0], 0, triangles_uv[i][3], triangles_uv[i][4], triangles_uv[i][5], zeta, fourier_coefficients, nfp, sine, cosine);
    //         z3 = coordinate_Z( 0, triangles_uv[i][2], triangles_uv[i][3], triangles_uv[i][4], triangles_uv[i][5], zeta, fourier_coefficients, nfp, sine, cosine);

    //         printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n",r1,z1,r2,z2,r3,z3);

    //         // temp1 += r1+z1+r2+z2+r3+z3;
    //     }
    // // }
    // // printf("temp1 = %f\n", temp1);
    // timing();
    // exit(0);




    double ***G = G_tri_20_18(triangles_uv);
    timing(time_s);
    double ***M1 = M_tri_18_18(G, triangles_uv);
    timing(time_s);
    // exit(0);
    double **GLs = Gauss_Legendre_quadrature_5();
    timing(time_s);
    // double t = numerical_integration_over_triangle(3.0, 5.0, 4.0, test2, GLs, 1);
    // printf("t = %0.32f, a = %0.32f\n",t, integral1(3, 1, 3.0, 5.0, 4.0));
    // exit(0);
    
    // int sections = num_sections;  // num_sections = 1 for test
    // double ***R_nu_sec_tri_18 = integral_functional1_nu_dudv(G, num_sections, GLs, (int) 1, coordinate_R, (int) nfp, fourier_coefficients, triangles_uv);
    // timing();
    // double ***Z_nu_sec_tri_18 = integral_functional1_nu_dudv(G, num_sections, GLs, (int) 1, coordinate_Z, (int) nfp, fourier_coefficients, triangles_uv);
    // timing();
    // double ***R_sec_tri_18 = F_sec_tri_18(M, R_nu_sec_tri_18, (int) num_sections, (int) triangles_uv[0][6]);
    // timing();
    // double ***Z_sec_tri_18 = F_sec_tri_18(M, Z_nu_sec_tri_18, (int) num_sections, (int) triangles_uv[0][6]);
    // timing();
    // //for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    // //    printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n", \
    // //    R_sec_tri_18[0][i][0],Z_sec_tri_18[0][i][0], R_sec_tri_18[0][i][6],Z_sec_tri_18[0][i][6], R_sec_tri_18[0][i][12],Z_sec_tri_18[0][i][12]);
    // //}
    // //exit(0);

    double ***pressure_nu_sec_tri_18 = integral_pressure_nu_dudv(G, num_sections, GLs, (int) GL_order, pressure_initialization, triangles_uv, initial_pressure);
    timing(time_s);
    double ***pressure_sec_tri_18 = F_sec_tri_18(M1, pressure_nu_sec_tri_18, (int) num_sections, (int) triangles_uv[0][6]);
    timing(time_s);
    for (int i = 0; i < num_sections; i++) {
         for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            delete [] pressure_nu_sec_tri_18[i][j];          
        }
        delete [] pressure_nu_sec_tri_18[i];
    }
    delete [] pressure_nu_sec_tri_18;    

    int num_variables = 2; // R  Z  Rzeta Zzeta  Rzetazeta Zzetazeta 
    double ****func_nu_sec_tri_18_num = integral_functional2_nu_dudv(G, num_sections, GLs, (int) GL_order, coordinate_funRZ, \
                                        (int) nfp, fourier_coefficients, fc_zeta, fc_zetazeta, (int) num_variables, triangles_uv);
    timing(time_s);
    double ****func_sec_tri_18_num = F_sec_tri_18_num(M1, func_nu_sec_tri_18_num, (int) num_sections, (int) triangles_uv[0][6], num_variables);
    timing(time_s);
    for (int i = 0; i < num_sections; i++) {
         for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            for (int k = 0; k < 18; k++) {
                delete [] func_nu_sec_tri_18_num[i][j][k];
            }
            delete [] func_nu_sec_tri_18_num[i][j];          
        }
        delete [] func_nu_sec_tri_18_num[i];
    }
    delete [] func_nu_sec_tri_18_num;
    
    double ***R_sec_tri_18 = new double **[num_sections];
    double ***Z_sec_tri_18 = new double **[num_sections];
    double ***Rzeta_sec_tri_18 = new double **[num_sections];
    double ***Zzeta_sec_tri_18 = new double **[num_sections];
    double ***Rzetazeta_sec_tri_18 = new double **[num_sections];
    double ***Zzetazeta_sec_tri_18 = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        R_sec_tri_18[i] = new double *[(int) triangles_uv[0][6]];
        Z_sec_tri_18[i] = new double *[(int) triangles_uv[0][6]];
        Rzeta_sec_tri_18[i] = new double *[(int) triangles_uv[0][6]];
        Zzeta_sec_tri_18[i] = new double *[(int) triangles_uv[0][6]];
        Rzetazeta_sec_tri_18[i] = new double *[(int) triangles_uv[0][6]];
        Zzetazeta_sec_tri_18[i] = new double *[(int) triangles_uv[0][6]];
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            R_sec_tri_18[i][j] = new double [18]();
            Z_sec_tri_18[i][j] = new double [18]();
            Rzeta_sec_tri_18[i][j] = new double [18]();
            Zzeta_sec_tri_18[i][j] = new double [18]();
            Rzetazeta_sec_tri_18[i][j] = new double [18]();
            Zzetazeta_sec_tri_18[i][j] = new double [18]();
            for (int k = 0; k < 18; k++) {
                R_sec_tri_18[i][j][k] = func_sec_tri_18_num[i][j][k][0];
                Z_sec_tri_18[i][j][k] = func_sec_tri_18_num[i][j][k][1];
                // Rzeta_sec_tri_18[i][j][k] = func_sec_tri_18_num[i][j][k][2];
                // Zzeta_sec_tri_18[i][j][k] = func_sec_tri_18_num[i][j][k][3];
                // Rzetazeta_sec_tri_18[i][j][k] = func_sec_tri_18_num[i][j][k][4];
                // Zzetazeta_sec_tri_18[i][j][k] = func_sec_tri_18_num[i][j][k][5];
            }
        }
    }

    double dzeta = 2 * pi / num_sections;
    double *zeta = new double [num_sections]();
    for (int i = 0; i < num_sections; i++) {
        zeta[i] = dzeta * i;
    }

    // dzeta d^2zeta;
    double **R_tri_18N = new double *[(int) triangles_uv[0][6]];
    double **Z_tri_18N = new double *[(int) triangles_uv[0][6]];
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        R_tri_18N[i] = new double [18*num_modes]();
        Z_tri_18N[i] = new double [18*num_modes]();
    }
    #pragma omp parallel for
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
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
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
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
    
    double ****divcoeff_sec_tri_3_var = divergence_free_constraints(triangles_uv, num_sections, num_modes, nfp_ft, R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18);
    timing(time_s);

    double ***J_nu_sec_tri_18 = integral_jacobian_nu_dudv(G, num_sections, GLs, (int) GL_order, jacobian, R_sec_tri_18, Z_sec_tri_18, triangles_uv);
    timing(time_s);
    double ***J_sec_tri_18 = F_sec_tri_18(M1, J_nu_sec_tri_18, (int) num_sections, (int) triangles_uv[0][6]);
    timing(time_s);
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            delete [] J_nu_sec_tri_18[i][j];
        }
        delete [] J_nu_sec_tri_18[i];
    }
    delete [] J_nu_sec_tri_18;
    
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < 18; j++) {
    //         printf("%0.9f ", J_sec_tri_18[0][i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n\n\n");
    // exit(0);


    // int rank = 9;
    // double *****P_nu_sec_tri_18_9_9 = integral_functional4_nu_dudv(G, num_sections, GLs, (int) 1, P_9_9, \
    // R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, Rzetazeta_sec_tri_18, Zzetazeta_sec_tri_18, triangles_uv, (int) rank);
    // timing();
    // double *****P_sec_tri_18_9_9 = F_sec_tri_18_m_m(M, P_nu_sec_tri_18_9_9, (int) num_sections, (int) triangles_uv[0][6], (int) rank);
    // timing();
    // for (int i = 0; i < num_sections; i++) {
    //      for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
    //         for (int k = 0; k < 18; k++) {
    //             for (int s = 0; s < 5; s++) {
    //                 delete [] P_nu_sec_tri_18_9_9[i][j][k][s];
    //             }
    //             delete [] P_nu_sec_tri_18_9_9[i][j][k];
    //         }
    //         delete [] P_nu_sec_tri_18_9_9[i][j];    
    //     }
    //     delete [] P_nu_sec_tri_18_9_9[i];
    // }
    // delete [] P_nu_sec_tri_18_9_9;    

    int num_scalar = 3;
    int num_operators = 11;
    double ******beltrami_nu_sec_tri_18_11_sca_sca = integral_beltrami_nu_dudv(G, num_sections, GLs, (int) GL_order, FuncBeltrami, \
    R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, Rzetazeta_sec_tri_18, Zzetazeta_sec_tri_18, triangles_uv, (int) num_scalar, beltrami_mu);
    timing(time_s);
    double ******beltrami_sec_tri_18_11_sca_sca = F_sec_tri_18_derivates_m_m(M1, beltrami_nu_sec_tri_18_11_sca_sca, (int) num_sections, (int) triangles_uv[0][6], (int) num_scalar, (int) num_operators);
    timing(time_s);
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
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

    // int num_poisson = 10;
    double ****poisson_nu_sec_tri_18_10 = integral_poisson_nu_dudv(G, num_sections, GLs, (int) GL_order, FuncPoisson, \
    R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, Rzetazeta_sec_tri_18, Zzetazeta_sec_tri_18, triangles_uv, (int) 10);
    timing(time_s);
    double ****poisson_sec_tri_18_10 = F_sec_tri_18_num(M1, poisson_nu_sec_tri_18_10, (int) num_sections, (int) triangles_uv[0][6], (int) 10);
    timing(time_s);
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            for (int k = 0; k < 18; k++) {
                delete [] poisson_nu_sec_tri_18_10[i][j][k];
            }
            delete [] poisson_nu_sec_tri_18_10[i][j];
        }
        delete [] poisson_nu_sec_tri_18_10[i];
    }
    delete [] poisson_nu_sec_tri_18_10;      

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



    
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < 18; j++) {
    //         for (int k = 0; k < 11; k++) {
    //             for (int s = 0; s < num_scalar; s++) {
    //                 for (int t = 0; t < num_scalar; t++) {
    //                     printf("%0.9f ", beltrami_sec_tri_18_11_sca_sca[0][i][j][k][s][t]);
    //                 }
    //             }
    //         }
    //         printf("\n");
    //     }
    //     printf("\n\n\n");
    // }
    // exit(0);

    // int num_beltrami = 11;
    // double ****beltrami_nu_sec_tri_18_11 = integral_functional5_nu_dudv(G, num_sections, GLs, (int) GL_order, FuncPoisson, \
    // R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, Rzetazeta_sec_tri_18, Zzetazeta_sec_tri_18, triangles_uv, (int) num_beltrami);
    // timing(time_s);
    // double ****beltrami_sec_tri_18_11 = F_sec_tri_18_num(M1, beltrami_nu_sec_tri_18_11, (int) num_sections, (int) triangles_uv[0][6], (int) num_beltrami);
    // timing(time_s);
    // for (int i = 0; i < num_sections; i++) {
    //      for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
    //         for (int k = 0; k < 18; k++) {
    //             delete [] beltrami_nu_sec_tri_18_11[i][j][k];
    //         }
    //         delete [] beltrami_nu_sec_tri_18_11[i][j];          
    //     }
    //     delete [] beltrami_nu_sec_tri_18_11[i];
    // }
    // delete [] beltrami_nu_sec_tri_18_11;



    // double temp_u = 0, temp_v = 0, temp_p = 0;
    // for (int s = 0; s < 5; s++) {
    //     for (int t = 0; t < 5; t++) {
    //         for (int i = 0; i < 20; i++) {
    //             for (int j = 0; j < 18; j++) {
    //                 temp_p += G[396][i][j] * power_nonnegative(temp_u, m_order[i]) * power_nonnegative(temp_v, n_order[i]) * P_sec_tri_18_5_5[0][396][j][s][t];
    //             }
    //         }
    //         printf("%0.9f, ", temp_p);
    //         temp_p = 0;
    //     }
    //     printf("\n");
    // }
    // // printf("\n");
    // // temp_u = -0.095093462;
    // // temp_v =  0.001615232;
    // // for (int s = 0; s < 5; s++) {
    // //     for (int t = 0; t < 5; t++) {
    // //         for (int i = 0; i < 20; i++) {
    // //             for (int j = 0; j < 18; j++) {
    // //                 temp_p += G[396][i][j] * power_nonnegative(temp_u, m_order[i]) * power_nonnegative(temp_v, n_order[i]) * P_sec_tri_18_9_9[0][396][j][s][t];
    // //             }
    // //         }
    // //         printf("%0.9f, ", temp_p);
    // //         temp_p = 0;
    // //     }
    // //     printf("\n");
    // // }
    // printf("\n");
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     printf("i = %d\n", i);
    //     for (int j = 0; j < 5; j++) {
    //         for (int k = 0; k < 5; k++) {
    //             for (int l = 0; l < 18; l++) {
    //                 printf("%0.9f, ", P_sec_tri_18_9_9[0][i][l][j][k]);
    //             }
    //             printf("\n");
    //         }
    //     }
    // }

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
    
    
    // double dzeta = 2 * pi / num_sections;
    // double *zeta = new double [num_sections]();
    // for (int i = 0; i < num_sections; i++) {
    //     zeta[i] = dzeta * i;
    // }

    double ***func_tri_18N_num = new double **[(int) triangles_uv[0][6]];
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        func_tri_18N_num[i] = new double *[18*num_modes];
        for (int j = 0; j < 18*num_modes; j++) {
            func_tri_18N_num[i][j] = new double [num_variables]();
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < num_variables; k++) {
                for (int l = 0; l < num_sections; l++) {
                    #pragma omp atomic
                    func_tri_18N_num[i][j][k] += func_sec_tri_18_num[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) *dzeta;
                }
            }
        }
    }
    for (int i = 0; i < num_sections; i++) {
         for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            for (int k = 0; k < 18; k++) {
                delete [] func_sec_tri_18_num[i][j][k];
            }
            delete [] func_sec_tri_18_num[i][j];          
        }
        delete [] func_sec_tri_18_num[i];
    }
    delete [] func_sec_tri_18_num;
    timing(time_s);  

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

    double **pressure_tri_18N = new double *[(int) triangles_uv[0][6]];
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        pressure_tri_18N[i] = new double [18*num_modes]();
        for (int j = 0; j < 18*num_modes; j++) {
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < num_sections; k++) {
                #pragma omp atomic
                pressure_tri_18N[i][j] += pressure_sec_tri_18[k][i][j/num_modes] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[k]) *dzeta;
            }
        }
    }
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            delete [] pressure_sec_tri_18[i][j];
        }
        delete [] pressure_sec_tri_18[i];
    }
    delete [] pressure_sec_tri_18;
    timing(time_s); 
    // print_p(pressure_tri_18N, G, triangles_uv, num_modes, nfp_ft);     
    // exit(0);

    double **J_tri_18N = new double *[(int) triangles_uv[0][6]];
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        J_tri_18N[i] = new double [18*num_modes]();
        for (int j = 0; j < 18*num_modes; j++) {
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < num_sections; k++) {
                #pragma omp atomic
                J_tri_18N[i][j] += J_sec_tri_18[k][i][j/num_modes] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[k]) *dzeta;
            }
        }
    }
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            delete [] J_sec_tri_18[i][j];
        }
        delete [] J_sec_tri_18[i];
    }
    delete [] J_sec_tri_18;
    timing(time_s);  


    double *****beltrami_tri_18N_11_sca_sca = new double ****[(int) triangles_uv[0][6]];
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
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
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < 11; k++) {
                for (int s = 0; s < num_scalar; s++) {
                    for (int t = 0; t < num_scalar; t++) {
                        for (int l = 0; l < num_sections; l++) {                            
                            #pragma omp atomic
                            beltrami_tri_18N_11_sca_sca[i][j][k][s][t] += beltrami_sec_tri_18_11_sca_sca[l][i][j/num_modes][k][s][t] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) *dzeta;
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
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


    double ***poisson_tri_18N_10 = new double **[(int) triangles_uv[0][6]];
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        poisson_tri_18N_10[i] = new double *[18*num_modes];
        for (int j = 0; j < 18*num_modes; j++) {
            poisson_tri_18N_10[i][j] = new double [10]();
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        for (int j = 0; j < 18*num_modes; j++) {
            for (int k = 0; k < 10; k++) {
                for (int l = 0; l < num_sections; l++) {                            
                    #pragma omp atomic
                    poisson_tri_18N_10[i][j][k] += poisson_sec_tri_18_10[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) *dzeta;
                }
            }
        }
    }
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
            for (int k = 0; k < 18; k++) {
                delete [] poisson_sec_tri_18_10[i][j][k];
            }
            delete [] poisson_sec_tri_18_10[i][j];
        }
        delete [] poisson_sec_tri_18_10[i];
    }
    delete [] poisson_sec_tri_18_10;    
    timing(time_s);    

    // double ****divcoeff_tri_3_var_modes = new double ***[(int) triangles_uv[0][6]];
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     divcoeff_tri_3_var_modes[i] = new double **[3];
    //     for (int j = 0; j < 3; j++) {
    //         divcoeff_tri_3_var_modes[i][j] = new double *[9];
    //         for (int k = 0; k < 9; k++) {
    //             divcoeff_tri_3_var_modes[i][j][k] = new double [num_modes]();
    //         }
    //     }
    // }
    // #pragma omp parallel for
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < 18*num_modes; j++) {
    //         for (int k = 0; k < 9; k++) {
    //             for (int l = 0; l < num_sections; l++) {                            
    //                 #pragma omp atomic
    //                 divcoeff_tri_3_var_modes[i][j/(6*num_modes)][k][j%num_modes] += divcoeff_sec_tri_18_9[l][i][j/num_modes][k] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) *dzeta;
    //             }
    //         }
    //     }
    // }
    // for (int i = 0; i < num_sections; i++) {
    //     for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
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
    double ***M2 = M_tri_18Nn_18Nn(MG, beltrami_tri_18N_11_sca_sca, M2_integral4, (int) triangles_uv[0][6], num_modes, nfp_ft, num_scalar, error, num_operators-1);
    timing(time_s);
    double ***M2_poisson = Mpoisson_tri_18N_18N(MG, poisson_tri_18N_10, M2_integral4, (int) triangles_uv[0][6], num_modes, nfp_ft, num_scalar, error);
    timing(time_s);
    // for (int i = 0; i < 1; i++) {
    //     for (int k = 0; k < 18*num_modes*num_scalar; k++) {
    //         for (int j = 0; j < 18*num_modes*num_scalar; j++) {
    //             printf("%0.9f ",M2_divergence[i][k][j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
    // }
    // exit(0);
    double **B2 = B_tri_18Nn(BG, beltrami_tri_18N_11_sca_sca, (int) triangles_uv[0][6], num_modes, nfp_ft, num_scalar, error, num_operators-1);
    timing(time_s);
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
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
    // solve_poisson_sparse(M2, B2, points_sequence, num_modes, (int) triangles_uv[0][6], time_s);
    // double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_direct(M2, B2, points_sequence, boundary_DoFs, num_modes, (int) triangles_uv[0][6], (int) num_scalar, error, time_s);
    // double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_lagrange_multiplier(M2, B2, points_sequence, boundary_DoFs, num_modes, (int) triangles_uv[0][6], (int) num_scalar, (int) 0, error, time_s);
    // double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_lagrange_multiplier_div(divcoeff_sec_tri_3_var, M2, B2, points_sequence, boundary_DoFs, num_modes, (int) triangles_uv[0][6], (int) nfp_ft, (int) num_scalar, (int) 0, error, time_s);
    // double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_lagrange_multiplier_div(divcoeff_sec_tri_3_var, M2, B2, points_sequence, normal_3_fixed, num_modes, (int) triangles_uv[0][6], (int) nfp_ft, (int) num_scalar, (int) 1, error, time_s);
    double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_ALM_div(divcoeff_sec_tri_3_var, M2, B2, points_sequence, boundary_DoFs, num_modes, (int) triangles_uv[0][6], (int) nfp_ft, (int) num_scalar, (int) 0, error, time_s);
    // double ***beltrami_3_tri_18N = solve_magnetic_field_sparse_ALM_div(divcoeff_sec_tri_3_var, M2, B2, points_sequence, normal_3_fixed, num_modes, (int) triangles_uv[0][6], (int) nfp_ft, (int) num_scalar, (int) 1, error, time_s);
    for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
        for (int j = 0; j < (18*num_modes)*num_scalar; j++) {
            delete [] M2[i][j];
        }
        delete [] M2[i];
        delete [] B2[i];
    }
    delete [] M2;
    delete [] B2;
    timing(time_s);
    // print_RZ(R_tri_18N, Z_tri_18N, G, triangles_uv, num_modes, nfp_ft);
    print_B(beltrami_3_tri_18N, G, triangles_uv, num_modes, nfp_ft, num_scalar);
    // exit(0);

    double ***mg1 = beltrami_3_tri_18N;
    double ***mg2;

    for (int i = 0; i < 2; i++) {
        mg2 = divergence_cleaning(mg1, M2_poisson, nfp_ft, num_modes, num_sections, \
        G, GLs, GL_order, R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, triangles_uv,  time_s, \
        M1, \
        MG, M2_integral4, BG, error, \
        num_segments, points_sequence, boundary_DoFs);
        print_B(mg2, G, triangles_uv, num_modes, nfp_ft, num_scalar);

        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < (int) triangles_uv[0][6]; k++) {
                for (int l = 0; l < 18*num_modes; l++) {
                    mg1[j][k][l] =  mg2[j][k][l];
                }
            }
        }

        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < (int) triangles_uv[0][6]; k++) {
                delete [] mg2[j][k];
            }
            delete [] mg2[j];
        }
        delete [] mg2;

        // mg2 = solve_magnetic_field_sparse_guess(mg1, M2, B2, points_sequence, boundary_DoFs, num_modes, (int) triangles_uv[0][6], (int) num_scalar, error, time_s);
        // print_B(mg2, G, triangles_uv, num_modes, nfp_ft, num_scalar);

        // for (int j = 0; j < 3; j++) {
        //     for (int k = 0; k < (int) triangles_uv[0][6]; k++) {
        //         for (int l = 0; l < 18*num_modes; l++) {
        //             mg1[j][k][l] =  mg2[j][k][l];
        //         }
        //     }
        // }

        // for (int j = 0; j < 3; j++) {
        //     for (int k = 0; k < (int) triangles_uv[0][6]; k++) {
        //         delete [] mg2[j][k];
        //     }
        //     delete [] mg2[j];
        // }
        // delete [] mg2;
    }
    // for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
    //     for (int j = 0; j < 18*num_modes; j++) {
    //         printf("%0.9f ", mg1[2][i][j]);
    //         if ((j+1) % (num_modes) == 0) {
    //             printf("         ");
    //         }
    //         if ((j+1) % (6*num_modes) == 0) {
    //             printf("\n");
    //         }
    //     }
    //     printf("\n");
    // }   
    exit(0);












    // SIME
    printf("solving for static ideal mhd equilibrium...\n");
    double ******mf_nu_sec_tri_18_5_sca_sca;
    double ******mf_sec_tri_18_5_sca_sca;
    double *****mf_tri_18N_5_sca_sca;
    for (int iter = 0; iter < 10; iter++) {
        mf_nu_sec_tri_18_5_sca_sca = integral_SMIEmf_nu_dudv(G, num_sections, GLs, (int) GL_order, Func_SIME_magnetic_field, \
        R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, mg1, pressure_tri_18N, triangles_uv, num_modes, nfp_ft, num_scalar);
        timing(time_s);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
                delete [] mg1[i][j];
            }
            delete [] mg1[i];
        }
        delete [] mg1;
        mf_sec_tri_18_5_sca_sca = F_sec_tri_18_derivates_m_m(M1, mf_nu_sec_tri_18_5_sca_sca, (int) num_sections, (int) triangles_uv[0][6], (int) num_scalar, int (5));
        timing(time_s);
        for (int i = 0; i < num_sections; i++) {
            for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
                for (int k = 0; k < 18; k++) {
                    for (int s = 0; s < 5; s++) {
                        for (int t = 0; t < num_scalar; t++) {
                            delete [] mf_nu_sec_tri_18_5_sca_sca[i][j][k][s][t];
                        }
                        delete [] mf_nu_sec_tri_18_5_sca_sca[i][j][k][s];
                    }
                    delete [] mf_nu_sec_tri_18_5_sca_sca[i][j][k];
                }
                delete [] mf_nu_sec_tri_18_5_sca_sca[i][j];
            }
            delete [] mf_nu_sec_tri_18_5_sca_sca[i];
        }
        delete [] mf_nu_sec_tri_18_5_sca_sca;  

        mf_tri_18N_5_sca_sca = new double ****[(int) triangles_uv[0][6]];
        for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
            mf_tri_18N_5_sca_sca[i] = new double ***[18*num_modes];
            for (int j = 0; j < 18*num_modes; j++) {
                mf_tri_18N_5_sca_sca[i][j] = new double **[5];
                for (int k = 0; k < 5; k++) {
                    mf_tri_18N_5_sca_sca[i][j][k] = new double *[num_scalar];
                    for (int s = 0; s < num_scalar; s++) {
                        mf_tri_18N_5_sca_sca[i][j][k][s] = new double [num_scalar]();
                    }
                }
            }
        }
        #pragma omp parallel for
        for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
            for (int j = 0; j < 18*num_modes; j++) {
                for (int k = 0; k < 5; k++) {
                    for (int s = 0; s < num_scalar; s++) {
                        for (int t = 0; t < num_scalar; t++) {
                            for (int l = 0; l < num_sections; l++) {                            
                                #pragma omp atomic
                                mf_tri_18N_5_sca_sca[i][j][k][s][t] += mf_sec_tri_18_5_sca_sca[l][i][j/num_modes][k][s][t] * fourier_basis(j%num_modes, (int) nfp_ft, zeta[l]) *dzeta;
                            }
                        }
                    }
                }
            }
        }
        for (int i = 0; i < num_sections; i++) {
            for (int j = 0; j < (int) triangles_uv[0][6]; j++) {
                for (int k = 0; k < 18; k++) {
                    for (int s = 0; s < 5; s++) {
                        for (int t = 0; t < num_scalar; t++) {
                            delete [] mf_sec_tri_18_5_sca_sca[i][j][k][s][t];
                        }
                        delete [] mf_sec_tri_18_5_sca_sca[i][j][k][s];
                    }
                    delete [] mf_sec_tri_18_5_sca_sca[i][j][k];
                }
                delete [] mf_sec_tri_18_5_sca_sca[i][j];
            }
            delete [] mf_sec_tri_18_5_sca_sca[i];
        }
        delete [] mf_sec_tri_18_5_sca_sca;    
        timing(time_s);   

        M2 = M_tri_18Nn_18Nn(MG, mf_tri_18N_5_sca_sca, M2_integral4, (int) triangles_uv[0][6], num_modes, nfp_ft, num_scalar, error, 4);
        timing(time_s);
        B2 = B_tri_18Nn(BG, mf_tri_18N_5_sca_sca, (int) triangles_uv[0][6], num_modes, nfp_ft, num_scalar, error, 4);
        timing(time_s);
        for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
            for (int j = 0; j < 18*num_modes; j++) {
                for (int k = 0; k < 5; k++) {
                    for (int l = 0; l < num_scalar; l++) {
                        delete [] mf_tri_18N_5_sca_sca[i][j][k][l];
                    }
                    delete [] mf_tri_18N_5_sca_sca[i][j][k];
                }
                delete [] mf_tri_18N_5_sca_sca[i][j];
            }
            delete [] mf_tri_18N_5_sca_sca[i];
        }
        delete [] mf_tri_18N_5_sca_sca;   

        // mg1 = solve_magnetic_field_sparse_direct(M2, B2, points_sequence, boundary_DoFs, num_modes, (int) triangles_uv[0][6], (int) num_scalar, error, time_s);
        mg1 = solve_magnetic_field_sparse_ALM_div(divcoeff_sec_tri_3_var, M2, B2, points_sequence, boundary_DoFs, num_modes, (int) triangles_uv[0][6], (int) nfp_ft, (int) num_scalar, (int) 0, error, time_s);
        for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
            for (int j = 0; j < (18*num_modes)*num_scalar; j++) {
                delete [] M2[i][j];
            }
            delete [] M2[i];
            delete [] B2[i];
        }
        delete [] M2;
        delete [] B2;
        timing(time_s);
        print_B(mg1, G, triangles_uv, num_modes, nfp_ft, num_scalar); 

        // for (int i = 0; i < 1; i++) {
        //     mg2 = divergence_cleaning(mg1, M2_poisson, nfp_ft, num_modes, num_sections, \
        //     G, GLs, GL_order, R_sec_tri_18, Z_sec_tri_18, Rzeta_sec_tri_18, Zzeta_sec_tri_18, triangles_uv,  time_s, \
        //     M1, \
        //     MG, M2_integral4, BG, error, \
        //     num_segments, points_sequence, boundary_DoFs);
        //     // print_B(mg2, G, triangles_uv, num_modes, nfp_ft, num_scalar);

        //     for (int j = 0; j < 2; j++) {
        //         for (int k = 0; k < (int) triangles_uv[0][6]; k++) {
        //             for (int l = 0; l < 18*num_modes; l++) {
        //                 mg1[j][k][l] =  mg2[j][k][l];
        //             }
        //         }
        //     }

        //     for (int j = 0; j < 3; j++) {
        //         for (int k = 0; k < (int) triangles_uv[0][6]; k++) {
        //             delete [] mg2[j][k];
        //         }
        //         delete [] mg2[j];
        //     }
        //     delete [] mg2;
        // }s


    }
    printf("\n\n\n\n\n----------------------------------------------\n");
    print_B(mg1, G, triangles_uv, num_modes, nfp_ft, num_scalar);
    




   
    printf("finished...\n");
    return 0;
}

