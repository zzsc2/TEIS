#include "infrastructure.hpp"

void timing() {
    time_t time;
    time = clock();
    printf(" time = %f s\n", double(time)/CLOCKS_PER_SEC);
    fflush(stdout);
}
void timing(time_t time_s) {
    time_t time_e;
    time_e = time(NULL);
    printf(" time = %d s\n", (int) (time_e - time_s));
    fflush(stdout);
}

template<class T>
T power_nonnegative(T x, int n) {
    if (n == 0) {
        return 1;
    }
    else if (n < 0) {
        return 0;
    }
    
    T value = 1;
    for (int i = 1; i <= n; i++) {
        value *= x;
    }
    return value;
}

double factorial(double n) {
    // return tgamma(n + 1.0);
    if ((int) n == 0) {
        return 1.0;
    }

    double value = 1;
    for (int i = 1; i <= (int) n; i++) {
        value *= i;
    }
    return value;
}

void table() {
    for (int i = 0; i < n_table; i++) {
        sine_table[i] = sin(step_r*i);
        cosine_table[i] = cos(step_r*i);
    }
}

double cosine_lookup_approximation(double angle) {
    // double x1 = fmod(abs(angle), 2*pi);
    double x1 = abs(angle);
    for (;;) {
        if (x1 < 2*pi) {
            break;
        }
        x1 -= 2*pi;
    }
    int i = (int) (x1 / step_r);
    double dx = x1 - (double) i * step_r;
    return cosine_table[i] * (1.0 - 0.5*dx*dx) - sine_table[i] * (dx);
}

double sine_lookup_approximation(double angle) {

    // double x1 = fmod(angle, 2*pi);
    // if (x1 < 0) {
    //     x1 += 2*pi;
    // }
    double x1 = angle;
    if (x1 >= 0) {
        for (;;) {
            if (x1 < 2*pi) {
                break;
            }
            x1 -= 2*pi;
        }
    }
    else {
        for (;;) {
            if (x1 >= 0) {
                break;
            }
            x1 += 2*pi;
        }
    }
    int i = (int) (x1 / step_r);
    double dx = x1 - (double) i * step_r;
    return sine_table[i] * (1.0 - 0.5*dx*dx) + cosine_table[i] * (dx);
}

double cosine_MMP_approximation(double angle) {
    // order = 16
    double x1 = fmod(abs(angle), 2*pi), x2, value;

    if (x1 >= pi) {
        x1 = 2*pi - x1;
        if (x1 >= pi/4) {
            x1 = pi - x1;
            x2 = x1 * x1;
            value = 0.99999999999999999608981951072301546 + x2*(-0.49999999999999974308667070690271678 + x2*(0.041666666666663887964497659102665977 + x2*(-0.00138888888887731722447063574018337351 + x2*(0.0000248015872774439536215582612189964674 + x2*(-2.75573163935355075808996593111486352e-7 + x2*(2.08765619601383980877940138085106396e-9 + x2*(-1.14629048996344469555493302470243688e-11 + 4.60900737685258733987906059894416831e-14*x2)))))));
            return -value;
        }
    }
    
    x2 = x1 * x1;
    value = 0.99999999999999999608981951072301546 + x2*(-0.49999999999999974308667070690271678 + x2*(0.041666666666663887964497659102665977 + x2*(-0.00138888888887731722447063574018337351 + x2*(0.0000248015872774439536215582612189964674 + x2*(-2.75573163935355075808996593111486352e-7 + x2*(2.08765619601383980877940138085106396e-9 + x2*(-1.14629048996344469555493302470243688e-11 + 4.60900737685258733987906059894416831e-14*x2)))))));
    return value;
}

double sine_MMP_approximation(double angle) {
    // order = 17
    double x1 = fmod(angle, 2*pi), x2, value;

    if (x1 < 0) {
        x1 += 2*pi;
    }

    if (x1 > 1.5*pi) {
        x1 = 2*pi - x1;
        x2 = x1 * x1;
        value = x1*(0.99999999999999999804111938371934874 + x2*(-0.166666666666666619002582182668768551 + x2*(0.00833333333333299311754500597897929703 + x2*(-0.000198412698411594497842939169564714446 + x2*(2.75573192045805570474164069181139307e-6 + x2*(-2.50521063810832659313336530795770274e-8 + x2*(1.60589186439698687926057420595617862e-10 + x2*(-7.64251167703190749594387242108845966e-13 + 2.71666540637857550793519752037048238e-15*x2))))))));
        return -value;
    }
    if (x1 > 1.0*pi && x1 <=1.5*pi) {
        x1 -= pi;
        x2 = x1 * x1;
        value = x1*(0.99999999999999999804111938371934874 + x2*(-0.166666666666666619002582182668768551 + x2*(0.00833333333333299311754500597897929703 + x2*(-0.000198412698411594497842939169564714446 + x2*(2.75573192045805570474164069181139307e-6 + x2*(-2.50521063810832659313336530795770274e-8 + x2*(1.60589186439698687926057420595617862e-10 + x2*(-7.64251167703190749594387242108845966e-13 + 2.71666540637857550793519752037048238e-15*x2))))))));
        return -value;
    }
    if (x1 > 0.5*pi && x1 <=1.0*pi) {
        x1 = pi - x1;
        x2 = x1 * x1;
        value = x1*(0.99999999999999999804111938371934874 + x2*(-0.166666666666666619002582182668768551 + x2*(0.00833333333333299311754500597897929703 + x2*(-0.000198412698411594497842939169564714446 + x2*(2.75573192045805570474164069181139307e-6 + x2*(-2.50521063810832659313336530795770274e-8 + x2*(1.60589186439698687926057420595617862e-10 + x2*(-7.64251167703190749594387242108845966e-13 + 2.71666540637857550793519752037048238e-15*x2))))))));
        return value;
    }
    else {
        x2 = x1 * x1;
        value = x1*(0.99999999999999999804111938371934874 + x2*(-0.166666666666666619002582182668768551 + x2*(0.00833333333333299311754500597897929703 + x2*(-0.000198412698411594497842939169564714446 + x2*(2.75573192045805570474164069181139307e-6 + x2*(-2.50521063810832659313336530795770274e-8 + x2*(1.60589186439698687926057420595617862e-10 + x2*(-7.64251167703190749594387242108845966e-13 + 2.71666540637857550793519752037048238e-15*x2))))))));
        return value;
    }
}

int binomial_coefficient(int n, int k) {
    if (k > n)
        return 0;
    if (k == 0 || k == n)
        return 1;
 
    return binomial_coefficient(n - 1, k - 1) + binomial_coefficient(n - 1, k);
}

template<class T>
bool value_comparison(T a, T b, double error) {
    if (abs(a-b) < error) {
        return true;
    }
    else {
        return false;
    }
}
template<class T>
T *enlarge_1d(T *old_array, int old_size, int new_size) {

    T *new_array = new T [new_size];

    for (int i = 0; i < old_size; i++) {
        new_array[i] = old_array[i];
    }
    delete [] old_array;

    return new_array;
}

template<class T>
T *reduce_1d(T *old_array, int old_size, int location, int reduced_size) {

    T *new_array = new T [old_size-reduced_size];

    for (int i = 0; i < old_size-reduced_size; i++) {
        if (i < location) {
            new_array[i] = old_array[i];
        }
        else if (i >= location && i < old_size-reduced_size) {
            new_array[i] = old_array[i+reduced_size];
        }
        else {
            printf("error in array reduction...\n");
            abort();
        }    
    }

    delete[] old_array;

    return new_array;
}

template<class T>
T **reduce_2d(T **old_array, int old_size, int location, int reduced_size, int dimension) {

    T **new_array = new T* [old_size-reduced_size];
    for (int i = 0; i < old_size-reduced_size; i++) {
        new_array[i] = new T [dimension];
    }

    for (int i = 0; i < old_size-reduced_size; i++) {
        if (i < location) {
            for (int j = 0; j < dimension; j++) {
                new_array[i][j] = old_array[i][j];
            }
            
        }
        else if (i >= location && i < old_size-reduced_size) {
            for (int j = 0; j < dimension; j++) {
                new_array[i][j] = old_array[i+reduced_size][j];
            }  
        }
        else {
            printf("error in array reduction...\n");
            abort();
        }    
    }

    for (int i = 0; i < old_size; i++) {
        delete [] old_array[i];
    }
    delete [] old_array;

    return new_array;
}

template<class T>
T **enlarge_2d(T **old_array, int old_size, int new_size, int dimension) {

    T **new_array = new T *[new_size];
    for (int i = 0; i < new_size; i++) {
        new_array[i] = new T [dimension];
    }

    for (int i = 0; i < old_size; i++) {
        for (int j = 0; j < dimension; j++) {
            new_array[i][j] = old_array[i][j];       
        }
        delete [] old_array[i];
    }
    delete [] old_array;
    
    return new_array;
}

template<class T>
T **insert_2d(T **old_array, int old_size, int location, int added_size, int dimension) {
  
    T **new_array = new T *[old_size + added_size];
    for (int i = 0; i < old_size + added_size; i++) {
        new_array[i] = new T [dimension];
    }

    for (int i = 0; i < old_size + added_size; i++) {
        for (int j = 0; j < dimension; j++) {
            
            if (i <= location) {
                new_array[i][j] = old_array[i][j];
            }
            else if (i > location && i <= location + added_size) {
                new_array[i][j] = 0;
            }
            else if (i > location + added_size && i < old_size + added_size) {
                new_array[i][j] = old_array[i-added_size][j];
            }
            else {
                printf("error in array insertion...\n");
                abort();
            }
        }
    }

    for (int i = 0; i < old_size; i++) {       
        delete [] old_array[i];
    }

    delete [] old_array;

    return new_array;
}

double fourier_basis(int i, int nfp, double zeta) {

    if (i < 0) {
        printf("i be must greater than 0 in FT...\n");
        exit(0);
    }
    else if (i == 0) {
        return 1 / (2 * pi);
    }
    else if (i % 2 == 1) {
        return cos(nfp * ((i+1)/2) * zeta) / pi;
    }
    else if (i % 2 == 0) {
        return sin(nfp * (i/2) * zeta) / pi;
    }
    else {
        printf("error in FT...\n");
        exit(0);
    }

}

double fourier_basis_test(int i, int nfp, double zeta) {

    if (i < 0) {
        printf("i must be greater than 0 in FT...\n");
        exit(0);
    }
    else if (i == 0) {
        return 1.0;
    }
    else if (i % 2 == 1) {
        return cos(nfp * ((i+1)/2) * zeta);
    }
    else if (i % 2 == 0) {
        return sin(nfp * (i/2) * zeta);
    }
    else {
        printf("error in FT...\n");
        exit(0);
    }

}

double fourier_basis_test2(int i, int nfp, double zeta) {

    if (i < 0) {
        printf("error in FT...\n");
        exit(0);
    }
    else if (i == 0) {
        return 0.0;
    }
    else if (i % 2 == 1) {
        return - (nfp * ((i+1)/2)) * sin(nfp * ((i+1)/2) * zeta);
    }
    else if (i % 2 == 0) {
        return (nfp * (i/2)) * cos(nfp * (i/2) * zeta);
    }
    else {
        printf("error in FT...\n");
        exit(0);
    }

}

double fourier_basis_test3(int i, int nfp, double zeta) {

    if (i < 0) {
        printf("error in FT...\n");
        exit(0);
    }
    else if (i == 0) {
        return 0.0;
    }
    else if (i % 2 == 1) {
        return - (nfp * ((i+1)/2)) * (nfp * ((i+1)/2)) * cos(nfp * ((i+1)/2) * zeta);
    }
    else if (i % 2 == 0) {
        return - (nfp * (i/2)) * (nfp * (i/2)) * sin(nfp * (i/2) * zeta);
    }
    else {
        printf("error in FT...\n");
        exit(0);
    }

}

void print_RZ(double **R, double **Z, double ***G, double **triangles_uv, int num_modes, int nfp_ft) {
    double rz1[2], rz2[2], rz3[2], uv1[20], uv2[20], uv3[20], r_modes[18], z_modes[18], zzeta;
    for (int l = 2; l < 3; l++) {
        zzeta = 0.0 + l * pi / 4.0;
        for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
            for (int j = 0; j < 18; j++) {
                r_modes[j] = 0;
                z_modes[j] = 0;
            }
            for (int i = 0; i < 2; i++) {
                rz1[i] = 0;
                rz2[i] = 0;
                rz3[i] = 0;
            }
            
            for (int j = 0; j < 18*num_modes; j++) {
                r_modes[j/num_modes] += R[i][j] * fourier_basis_test(j%num_modes, nfp_ft, zzeta);
                z_modes[j/num_modes] += Z[i][j] * fourier_basis_test(j%num_modes, nfp_ft, zzeta);
            }
            for (int j = 0; j < 20; j++) {
                uv1[j] = power_nonnegative(-triangles_uv[i][1], m_order[j]) *   power_nonnegative(0, n_order[j]);
                uv2[j] = power_nonnegative( triangles_uv[i][0], m_order[j]) *   power_nonnegative(0, n_order[j]);
                uv3[j] = power_nonnegative(0, m_order[j]) *   power_nonnegative(triangles_uv[i][2], n_order[j]);
                for (int k = 0; k < 18; k++) {
                    rz1[0] += G[i][j][k] * uv1[j] * r_modes[k];
                    rz2[0] += G[i][j][k] * uv2[j] * r_modes[k];
                    rz3[0] += G[i][j][k] * uv3[j] * r_modes[k];

                    rz1[1] += G[i][j][k] * uv1[j] * z_modes[k];
                    rz2[1] += G[i][j][k] * uv2[j] * z_modes[k];
                    rz3[1] += G[i][j][k] * uv3[j] * z_modes[k];
                }
            }
            printf("%0.9f %0.9f ; %0.9f %0.9f ; %0.9f %0.9f;\n", rz1[0],rz1[1],rz2[0],rz2[1],rz3[0],rz3[1]);
        }
        printf("\n\n\n\n\n");
    }
    fflush(stdout);
}

void print_B(double ***F, double ***G, double **triangles_uv, int num_modes, int nfp_ft, int num_scalar) {

    double f1[num_scalar], f2[num_scalar], f3[num_scalar], uv1[20], uv2[20], uv3[20], f_modes_2d[num_scalar][18], zzeta;
    for (int l = 2; l < 3; l++) {
        zzeta = 0.0 + l * pi / 4.0;
        for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
            for (int j = 0; j < 18; j++) {
                for (int k = 0; k < num_scalar; k++) {
                    f_modes_2d[k][j] = 0;
                }
            }
            for (int j = 0; j < num_scalar; j++) {
                f1[j] = 0;
                f2[j] = 0;
                f3[j] = 0;
            }
            
            for (int k = 0; k < num_scalar; k++) {
                for (int j = 0; j < 18*num_modes; j++) {
                    f_modes_2d[k][j/num_modes] += F[k][i][j] * fourier_basis_test(j%num_modes, nfp_ft, zzeta);
                }
            }
            for (int j = 0; j < 20; j++) {
                uv1[j] = power_nonnegative(-triangles_uv[i][1], m_order[j]) *   power_nonnegative(0, n_order[j]);
                uv2[j] = power_nonnegative( triangles_uv[i][0], m_order[j]) *   power_nonnegative(0, n_order[j]);
                uv3[j] = power_nonnegative(0, m_order[j]) *   power_nonnegative(triangles_uv[i][2], n_order[j]);
                for (int k = 0; k < 18; k++) {
                    for (int s = 0; s < num_scalar; s++) {
                        f1[s] += G[i][j][k] * uv1[j] * f_modes_2d[s][k];
                        f2[s] += G[i][j][k] * uv2[j] * f_modes_2d[s][k];
                        f3[s] += G[i][j][k] * uv3[j] * f_modes_2d[s][k];
                    }
                }
            }
            if (num_scalar == 3) {
                printf("%0.9f  %0.9f  %0.9f\n", sqrt(f1[0]*f1[0]+f1[1]*f1[1]+f1[2]*f1[2]), \
                sqrt(f2[0]*f2[0]+f2[1]*f2[1]+f2[2]*f2[2]), \
                sqrt(f3[0]*f3[0]+f3[1]*f3[1]+f3[2]*f3[2]));
                // printf("%0.9f  %0.9f  %0.9f;    %0.9f  %0.9f  %0.9f;    %0.9f  %0.9f  %0.9f\n", \
                // f1[0], f2[0], f3[0], f1[1], f2[1], f3[1], f1[2], f2[2], f3[2]);
            }
            else if (num_scalar == 1) {
                printf("%0.9f  %0.9f  %0.9f\n", abs(f1[0]), abs(f2[0]), abs(f3[0]));
            }
        }
        printf("\n\n\n\n\n");
    }
    fflush(stdout);
}

void print_p(double **F, double ***G, double **triangles_uv, int num_modes, int nfp_ft) {

    double f1, f2, f3, uv1[20], uv2[20], uv3[20], f_modes[18], zzeta;
    for (int l = 2; l < 3; l++) {
        zzeta = 0.0 + l * pi / 4.0;
        for (int i = 0; i < (int) triangles_uv[0][6]; i++) {
            for (int j = 0; j < 18; j++) {
                f_modes[j] = 0;
            }
            f1 = 0;
            f2 = 0;
            f3 = 0;
            
            for (int j = 0; j < 18*num_modes; j++) {
                f_modes[j/num_modes] += F[i][j] * fourier_basis_test(j%num_modes, nfp_ft, zzeta);
            }
            for (int j = 0; j < 20; j++) {
                uv1[j] = power_nonnegative(-triangles_uv[i][1], m_order[j]) *   power_nonnegative(0, n_order[j]);
                uv2[j] = power_nonnegative( triangles_uv[i][0], m_order[j]) *   power_nonnegative(0, n_order[j]);
                uv3[j] = power_nonnegative(0, m_order[j]) *   power_nonnegative(triangles_uv[i][2], n_order[j]);
                for (int k = 0; k < 18; k++) {
                    f1 += G[i][j][k] * uv1[j] * f_modes[k];
                    f2 += G[i][j][k] * uv2[j] * f_modes[k];
                    f3 += G[i][j][k] * uv3[j] * f_modes[k];
                }
            }
            printf("%0.9f  %0.9f  %0.9f\n", f1, f2, f3);
        }
        printf("\n\n\n\n\n");
    }
    fflush(stdout);
}

double low_of_cosines(double a, double b, double c) {
    return acos((a*a+b*b-c*c) / (2*b*c));
}
