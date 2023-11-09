#include "integrals.hpp"

void mn_order() {
    m_order[0] = 0;
    m_order[1] = 1;
    m_order[2] = 0;
    m_order[3] = 2;
    m_order[4] = 1;
    m_order[5] = 0;
    m_order[6] = 3;
    m_order[7] = 2;
    m_order[8] = 1;
    m_order[9] = 0;
    m_order[10] = 4;
    m_order[11] = 3;
    m_order[12] = 2;
    m_order[13] = 1;
    m_order[14] = 0;
    m_order[15] = 5;
    m_order[16] = 3;
    m_order[17] = 2;
    m_order[18] = 1;
    m_order[19] = 0;

    n_order[0] = 0;
    n_order[1] = 0;
    n_order[2] = 1;
    n_order[3] = 0;
    n_order[4] = 1;
    n_order[5] = 2;
    n_order[6] = 0;
    n_order[7] = 1;
    n_order[8] = 2;
    n_order[9] = 3;
    n_order[10] = 0;
    n_order[11] = 1;
    n_order[12] = 2;
    n_order[13] = 3;
    n_order[14] = 4;
    n_order[15] = 0;
    n_order[16] = 2;
    n_order[17] = 3;
    n_order[18] = 4;
    n_order[19] = 5;
}

double integral1(int M, int N, double a, double b, double c) {
    if (M < 0 || N < 0) {
        return 0;
    }
    double m = (double) M;
    double n = (double) N;
    return factorial(m) * factorial(n) / factorial(m+n+2) * (power_nonnegative(a,M+1) - power_nonnegative(-b,M+1)) * power_nonnegative(c, N+1);
}

double integral2(double a, double b, double c, double ***G, int i, int j, int position) {
    double m_ij = 0;
    for (int k = 0; k < 20; k++) {
        for (int l = 0; l < 20; l++) {
            m_ij += G[position][k][i] * G[position][l][j] * integral1(m_order[k] + m_order[l], n_order[k] + n_order[l], a, b, c);
        }
    }
    return m_ij;
}


double integral3(double ***G, int i, int j, int k, int position, int m, int n, double *****array_integral3) {
    double m_ijk = 0;
    double temp = 0;
    for (int l = 0; l < 20; l++) {
        for (int s = 0; s < 20; s++) {
            for (int t = 0; t < 20; t++) {

                temp = G[position][l][i] * G[position][s][j] * G[position][t][k] * array_integral3[l][s][t][-m][-n];

                for (int i = 0; i < -m; i++) {
                    temp *= (m_order[t] - i);
                }
                for (int i = 0; i < -n; i++) {
                    temp *= (n_order[t] - i);
                }
                m_ijk += temp;

            }
        }
    }
    return m_ijk;
}

double integral4(int nfp, int ntor, int t, int k, int j, int Gj, double error) { 
    // G*G*G

    std::function<double(double, double, double)> sign_max = [&] (double x, double y, double z) {
        if (x > y && x > z) {
            return -1;
        }
        else {
            return 1;
        }
    };

    std::function<double(int, int, int, int, int, double)> func1 = [&](int nfp, int ntor, int t, int k, int j, double error) {
        // Gt = 0, Gk = 0, Gj = 0

        double value = 0;

        int k1 = t % ntor; 
        int k2 = k % ntor; 
        int k3 = j % ntor;

        int num1 = k1 % 2;
        int num2 = k2 % 2;
        int num3 = k3 % 2;

        double k1_0 = nfp * k1 / 2.0, k1_1 = nfp * (k1+1) / 2.0;
        double k2_0 = nfp * k2 / 2.0, k2_1 = nfp * (k2+1) / 2.0;
        double k3_0 = nfp * k3 / 2.0, k3_1 = nfp * (k3+1) / 2.0;

        double a, b, c;

        if (k1 != 0 && k2 != 0 && k3 != 0) {
            if (num1 == 0 && num2 == 0 && num3 == 1) {
                a = k1_0;
                b = k2_0;
                c = k3_1;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = sign_max(k3_1, k2_0, k1_0) * pi / 2.0;
                }
            }
            else if (num1 == 0 && num2 == 1 && num3 == 0) {
                a = k1_0;
                b = k2_1;
                c = k3_0;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = sign_max(k2_1, k3_0, k1_0) * pi / 2.0;
                }
            }
            else if (num1 == 1 && num2 == 0 && num3 == 0) {
                a = k1_1;
                b = k2_0;
                c = k3_0;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = sign_max(k1_1, k2_0, k3_0) * pi / 2.0;
                }
            }
            else if (num1 == 1 && num2 == 1 && num3 == 1) {
                a = k1_1;
                b = k2_1;
                c = k3_1;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = pi / 2.0;
                }
            }
        }
        else if (k1 == 0 && k2 != 0 && k3 != 0) {
            if (num2 == 0 && num3 == 0) {
                b = k2_0;
                c = k3_0;
                if (value_comparison((b-c), 0.0, error)) {
                    value = pi;
                }
            }
            else if (num2 == 1 && num3 == 1) {
                b = k2_1;
                c = k3_1;
                if (value_comparison((b-c), 0.0, error)) {
                    value = pi;
                }
            }
        }
        else if (k1 != 0 && k2 == 0 && k3 != 0) {
            if (num1 == 0 && num3 == 0) {
                b = k1_0;
                c = k3_0;
                if (value_comparison((b-c), 0.0, error)) {
                    value = pi;
                }
            }
            else if (num1 == 1 && num3 == 1) {
                b = k1_1;
                c = k3_1;
                if (value_comparison((b-c), 0.0, error)) {
                    value = pi;
                }
            }
        }
        else if (k1 != 0 && k2 != 0 && k3 == 0) {
            if (num1 == 0 && num2 == 0) {
                b = k1_0;
                c = k2_0;
                if (value_comparison((b-c), 0.0, error)) {
                    value = pi;
                }
            }
            else if (num1 == 1 && num2 == 1) {
                b = k1_1;
                c = k2_1;
                if (value_comparison((b-c), 0.0, error)) {
                    value = pi;
                }
            }
        }
        else if (k1 == 0 && k2 == 0 && k3 == 0) {
            value = 2.0 * pi;
        }
           
        return value;
    };

    std::function<double(int, int, int, int, int, double)> func2 = [&](int nfp, int ntor, int t, int k, int j, double error) {
        // Gt = 0, Gk = 0, Gj = 1

        double value = 0;

        int k1 = t % ntor; 
        int k2 = k % ntor; 
        int k3 = j % ntor;

        int num1 = k1 % 2;
        int num2 = k2 % 2;
        int num3 = k3 % 2;

        double k1_0 = nfp * k1 / 2.0, k1_1 = nfp * (k1+1) / 2.0;
        double k2_0 = nfp * k2 / 2.0, k2_1 = nfp * (k2+1) / 2.0;
        double k3_0 = nfp * k3 / 2.0, k3_1 = nfp * (k3+1) / 2.0;

        double a, b, c;

        if (k1 != 0 && k2 != 0 && k3 != 0) {
            if (num1 == 0 && num2 == 0 && num3 == 0) {
                a = k1_0;
                b = k2_0;
                c = k3_0;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = sign_max(k3_0, k2_0, k1_0) * k3_0 * pi / 2.0;
                }
            }
            else if (num1 == 0 && num2 == 1 && num3 == 1) {
                a = k1_0;
                b = k2_1;
                c = k3_1;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = - sign_max(k2_1, k3_1, k1_0) * k3_1 * pi / 2.0;
                }
            }
            else if (num1 == 1 && num2 == 0 && num3 == 1) {
                a = k1_1;
                b = k2_0;
                c = k3_1;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = - sign_max(k1_1, k2_0, k3_1) * k3_1 * pi / 2.0;
                }
            }
            else if (num1 == 1 && num2 == 1 && num3 == 0) {
                a = k1_1;
                b = k2_1;
                c = k3_0;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = k3_0 * pi / 2.0;
                }
            }
        }
        else if (k1 == 0 && k2 != 0 && k3 != 0) {
            if (num2 == 0 && num3 == 1) {
                b = k2_0;
                c = k3_1;
                if (value_comparison((b-c), 0.0, error)) {
                    value = - k3_1 * pi;
                }
            }
            else if (num2 == 1 && num3 == 0) {
                b = k2_1;
                c = k3_0;
                if (value_comparison((b-c), 0.0, error)) {
                    value = k3_0 * pi;
                }
            }
        }
        else if (k1 != 0 && k2 == 0 && k3 != 0) {
            if (num1 == 0 && num3 == 1) {
                b = k1_0;
                c = k3_1;
                if (value_comparison((b-c), 0.0, error)) {
                    value = - k3_1 * pi;
                }
            }
            else if (num1 == 1 && num3 == 0) {
                b = k1_1;
                c = k3_0;
                if (value_comparison((b-c), 0.0, error)) {
                    value = k3_0 * pi;
                }
            }
        }
          
        return value;
    };

    std::function<double(int, int, int, int, int, double)> func3 = [&](int nfp, int ntor, int t, int k, int j, double error) {
        // Gt = 0, Gk = 0, Gj = 2

        double value = 0;

        int k1 = t % ntor; 
        int k2 = k % ntor; 
        int k3 = j % ntor;

        int num1 = k1 % 2;
        int num2 = k2 % 2;
        int num3 = k3 % 2;

        double k1_0 = nfp * k1 / 2.0, k1_1 = nfp * (k1+1) / 2.0;
        double k2_0 = nfp * k2 / 2.0, k2_1 = nfp * (k2+1) / 2.0;
        double k3_0 = nfp * k3 / 2.0, k3_1 = nfp * (k3+1) / 2.0;

        double a, b, c;
        
        if (k1 != 0 && k2 != 0 && k3 != 0) {
            if (num1 == 0 && num2 == 0 && num3 == 1) {
                a = k1_0;
                b = k2_0;
                c = k3_1;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = - sign_max(k3_1, k2_0, k1_0) * k3_1 * k3_1 * pi / 2.0;
                }
            }
            else if (num1 == 0 && num2 == 1 && num3 == 0) {
                a = k1_0;
                b = k2_1;
                c = k3_0;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = - sign_max(k2_1, k3_0, k1_0) * k3_0 * k3_0 * pi / 2.0;
                }
            }
            else if (num1 == 1 && num2 == 0 && num3 == 0) {
                a = k1_1;
                b = k2_0;
                c = k3_0;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = - sign_max(k1_1, k2_0, k3_0) * k3_0 * k3_0 * pi / 2.0;
                }
            }
            else if (num1 == 1 && num2 == 1 && num3 == 1) {
                a = k1_1;
                b = k2_1;
                c = k3_1;
                if (value_comparison((a+b-c)*(c+a-b)*(b+c-a), 0.0, error)) {
                    value = - k3_1 * k3_1 * pi / 2.0;
                }
            }
        }
        else if (k1 == 0 && k2 != 0 && k3 != 0) {
            if (num2 == 0 && num3 == 0) {
                b = k2_0;
                c = k3_0;
                if (value_comparison((b-c), 0.0, error)) {
                    value = - k3_0 * k3_0 * pi;
                }
            }
            else if (num2 == 1 && num3 == 1) {
                b = k2_1;
                c = k3_1;
                if (value_comparison((b-c), 0.0, error)) {
                    value = - k3_1 * k3_1 * pi;
                }
            }
        }
        else if (k1 != 0 && k2 == 0 && k3 != 0) {
            if (num1 == 0 && num3 == 0) {
                b = k1_0;
                c = k3_0;
                if (value_comparison((b-c), 0.0, error)) {
                    value = - k3_0 * k3_0 * pi;
                }
            }
            else if (num1 == 1 && num3 == 1) {
                b = k1_1;
                c = k3_1;
                if (value_comparison((b-c), 0.0, error)) {
                    value = - k3_1 * k3_1 * pi;
                }
            }
        }
        
        return value;
    };

    if (Gj == 0) {
        return func1(nfp, ntor, t, k, j, error);
    }
    else if (Gj == 1) {
        return func2(nfp, ntor, t, k, j, error);
    }
    else if (Gj == 2) {
        return func3(nfp, ntor, t, k, j, error);
    }
    else {
        printf("integral3: wrong order of Gj...\n");
        exit(0);
    }
    
}


double **Gauss_Legendre_quadrature_5() {
    printf("calculating Gauss_Legendre_coefficients...\n");
    double GL_p[5] = {- sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0, \
                      - sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0, \
                        0, \
                        sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0, \
                        sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0};
    
    double GL_w[5] = {(322.0 - 13.0 * sqrt(70.0)) / 900.0, \
                      (322.0 + 13.0 * sqrt(70.0)) / 900.0, \
                       128.0 / 225.0, \
                      (322.0 + 13.0 * sqrt(70.0)) / 900.0, \
                      (322.0 - 13.0 * sqrt(70.0)) / 900.0};
    
    double c[25], x[25], y[25];
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            c[5*i+j] = (2.0 - GL_p[i] - GL_p[j]) * GL_w[i] * GL_w[j] / 4.0;
            x[5*i+j] = (- 1.0 + 3.0 * GL_p[i] - GL_p[j] * (1 + GL_p[i])) / 4.0;
            y[5*i+j] = (- 1.0 + 3.0 * GL_p[j] - GL_p[i] * (1 + GL_p[j])) / 4.0;
        }
    }

    double ** GLs = new double *[3];
    for (int i = 0; i < 3; i++) {
        GLs[i] = new double [25];
    }
    for (int i = 0; i < 25; i++) {
        GLs[0][i] = c[i];
        GLs[1][i] = x[i];
        GLs[2][i] = y[i];
        // printf("%0.9f %0.9f %0.9f\n", GLs[0][i], GLs[1][i], GLs[2][i]);
    }
    // exit(0);

    return GLs;
}

double numerical_integration_over_triangle(double A, double B, double C, std::function<double(double, double)> func, double **GLs, int n) {
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

    double h[25], val = 0;
    for (int k = 0; k < 25; k++) {
        h[k] = 0;
        for (int i = 0; i < 2 * n; i++) {   
            for (int j = 0; j < 2 * n - i; j++) {
                uu = (x[k] + 2*(i-n) + 1) / (2*n);
                vv = (y[k] + 2*(j-n) + 1) / (2*n);
                l1 = -0.5 * (uu + vv);
                l2 = 0.5 * (1 + uu);
                l3 = 0.5 * (1 + vv);
                u = -B * l1 + A *l2;
                v = C * l3;
                h[k] += func(u, v);
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
                h[k] += func(u, v);
            }
        }
        val += c[k] * h[k];  
    }
    val /= (4*n*n);
    val *= (A+B)*C/4;

    return val;
}

double numerical_integration1_over_triangle(double A, double B, double C, std::function<double(double, double, double, double, double, double, double**, int, double, double, int, int)> func, \
double theta, double x0, double y0, double zeta, double **fourier_coefficients, int nfp, double sine, double cosine, int position, int column, \
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
                h[k] += func(u, v, theta, x0, y0, zeta, fourier_coefficients, nfp, sine, cosine, position, column);
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
                h[k] += func(u, v, theta, x0, y0, zeta, fourier_coefficients, nfp, sine, cosine, position, column);
            }
        }

        val += c[k] * h[k];  

    }

    val /= (4*n*n);
    val *= (A+B)*C/4;

    return val;
}

double **numerical_integration4_over_triangle(double A, double B, double C, std::function<double**(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int, int, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri, int column, int rank, \
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

    double ***h = new double **[25];
    for (int i = 0; i < 25; i++) {
        h[i] = new double *[rank];
        for (int j = 0; j < rank; j++) {
            h[i][j] = new double [rank];
            for (int k = 0; k < rank; k++) {
                h[i][j][k] = 0;
            }
        }
    }

    double **val = new double *[rank];
    for (int i = 0; i < rank; i++) {
        val[i] = new double [rank];
        for (int j = 0; j < rank; j++) {
            val[i][j] = 0;
        }
    }

    double **ptr;
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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, position_sec, position_tri, column, rank);
                for (int s = 0; s < rank; s++) {
                    for (int t = 0; t < rank; t++) {
                        h[k][s][t] += ptr[s][t];
                    }
                }
                for (int l = 0; l < rank; l++) {
                    delete[] ptr[l];
                }
                delete[] ptr;
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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, position_sec, position_tri, column, rank);
                for (int s = 0; s < rank; s++) {
                    for (int t = 0; t < rank; t++) {
                        h[k][s][t] += ptr[s][t];
                    }
                }
                for (int l = 0; l < rank; l++) {
                    delete[] ptr[l];
                }
                delete[] ptr;
            }
        }

        for (int s = 0; s < rank; s++) {
            for (int t = 0; t < rank; t++) {
                val[s][t] += c[k] * h[k][s][t] ;
            }
        }
          
    }

    for (int s = 0; s < rank; s++) {
        for (int t = 0; t < rank; t++) {
            val[s][t] /= (4*n*n);
            val[s][t] *= (A+B)*C/4;
        }
    }

    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < rank; j++) {
            delete [] h[i][j];
        }
        delete [] h[i];
    }
    delete [] h;
// printf("position_sec = %d; position_tri = %d; column = %d\n",position_sec, position_tri, column);
    return val;
}



double ***integral_functional1_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double(double, double, double, double, double, double, double **, int, double, double)> functional, \
int nfp, double **fourier_coefficients, double **triangles_uv) {

    printf("calculating integral of functional1*nu_i...\n");
    int num_triangles_in_plane = (int) triangles_uv[0][6];
    
    std::function<double(double, double, double, double, double, double, double**, int, double, double, int, int)> functional_nu = \
    [&](double u, double v, double theta, double x0, double y0, double zeta, double **fc, int nfp, double sine, double cosine, int position, int column) {

        double nu_i = 0, value;
        for (int j = 0; j < 20; j++) {
            nu_i += (G[position][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }
        value = functional(u, v, theta, x0, y0, zeta, fc, nfp, sine, cosine) * nu_i;

        return value;
    };  

    double ***functional_nu_sec_tri_18 = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18[i] = new double *[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            // printf("section: %d; triangle: %d\n", i, j);
            functional_nu_sec_tri_18[i][j] = new double [18];            
        }         
    }

    for (int i = 0; i < num_sections; i++) {
        #pragma omp parallel for num_threads (NT)
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {
                // column = k;
                // functional_nu_sec_tri_18[i][j][k] = numerical_integration_over_triangle(a, b, c, functional_nu, \
                // theta, x0, y0, zeta, fourier_coefficients, nfp, sine, cosine, position, column, \
                // GLs, n);

                functional_nu_sec_tri_18[i][j][k] = numerical_integration1_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                triangles_uv[j][3], triangles_uv[j][4], triangles_uv[j][5], (double) i * (2*pi / num_sections), fourier_coefficients, nfp, sin(triangles_uv[j][3]), cos(triangles_uv[j][3]), j, k, \
                GLs, n);
            } 
        }
    }

    return functional_nu_sec_tri_18;
}

double *****integral_functional4_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double **(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double **triangles_uv, int rank) {

    printf("calculating integral of P*nu_i...\n");
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double**(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int, int, int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri, int column, int rank) {

        double nu_i = 0;

        double **value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, position_sec, position_tri);
        for (int i = 0; i < rank; i++) {
            for (int j = 0; j < rank; j++) {
                value[i][j] *= nu_i;
            }
        }

        return value;
    };

    double *****functional_nu_sec_tri_18_m_m = new double ****[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_m_m[i] = new double ***[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_m_m[i][j] = new double **[18];
        }
    }

    for (int i = 0; i < num_sections; i++) {
        #pragma omp parallel for num_threads (NT)
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_m_m[i][j][k] = numerical_integration4_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, i, j, k, rank, \
                                                                                             GLs, n);
            } 
        }
    }

    return functional_nu_sec_tri_18_m_m;
}


//
double *numerical_integration_DivCoeff_over_triangle(double A, double B, double C, std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int column, \
double **GLs, int n) {
    // order = 5

    int num = 9;

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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri, column);
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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri, column);
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

double ****integral_DivCoeff_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double **triangles_uv) {

    printf("calculating integral of (DivCoeff)*nu_i...\n");
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int column) {

        double nu_i = 0;

        int num_variables = 9;

        double *value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri);
        for (int i = 0; i < num_variables; i++) {
            value[i] *= nu_i;
        }

        return value;
    };

    double ****functional_nu_sec_tri_18_m = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_m[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_m[i][j] = new double *[18];
        }
    }

    for (int i = 0; i < num_sections; i++) {
        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_m[i][j][k] = numerical_integration_DivCoeff_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, G, i, j, k, \
                                                                                             GLs, n);
            } 
        }
    }

    return functional_nu_sec_tri_18_m;
}

//
double *numerical_integration_projection_over_triangle(double A, double B, double C, std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int, int, int, double ***, double ***, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ***B_tri_3_18, double ***poisoon_2_tri_18, int column, \
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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri, num_modes, nfp_ft, num_sections, B_tri_3_18, poisoon_2_tri_18, column);
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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri, num_modes, nfp_ft, num_sections, B_tri_3_18, poisoon_2_tri_18, column);
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

double ****integral_projection_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int, int, int, double ***, double ***)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double **triangles_uv, int num_modes, int nfp_ft, double ***B_3_tri_18N, double ***poisoon_1_tri_18N) {

    printf("calculating integral of (Projection)*nu_i...\n");
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, int, int, int, int, int, double ***, double ***, int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ***B_tri_3_18, double ***poisoon_2_tri_18, int column) {

        double nu_i = 0;

        int num = 3;

        double *value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri, num_modes, nfp_ft, num_sections, B_tri_3_18, poisoon_2_tri_18);
        for (int i = 0; i < num; i++) {
            value[i] *= nu_i;
        }

        return value;
    };

    double ****functional_nu_sec_tri_18_m = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_m[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_m[i][j] = new double *[18];
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

        // for (int i = 0; i < num_triangles_in_plane; i++) {
        //     for (int j = 0; j < 18; j++) {
        //         printf("%0.9f ", poisoon_2_tri_18[1][i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n\n\n\n");
        // fflush(stdout);
        // exit(0);

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_m[i][j][k] = numerical_integration_projection_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
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

    return functional_nu_sec_tri_18_m;
}

//
double numerical_integration_divergence_over_triangle(double A, double B, double C, std::function<double(double , double , double ***, double ***, double ***, double ***, double ***, int , int , int , int , int, double ****, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ****B_2_tri_3_18, int column, \
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

    double h[25] = {0};


    double val = 0.0;

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
                h[k] += func(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri, num_modes, nfp_ft, num_sections, B_2_tri_3_18, column);
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
                h[k] += func(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri, num_modes, nfp_ft, num_sections, B_2_tri_3_18, column);
            }
        }

        val += c[k] * h[k];      

    }
   
    val /= (4*n*n);
    val *= (A+B)*C/4;

    return val;
}

double ***integral_divergence_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double (double , double , double ***, double ***, double ***, double ***, double ***, int , int , int , int , int, double ****)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double **triangles_uv, int num_modes, int nfp_ft, double ***B_3_tri_18N) {

    printf("calculating integral of (Divergence)*nu_i...\n");
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double(double , double , double ***, double ***, double ***, double ***, double ***, int , int , int , int , int, double ****, int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ****B_2_tri_3_18, int column) {

        double nu_i = 0;

        double value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, G, position_sec, position_tri, num_modes, nfp_ft, num_sections, B_2_tri_3_18) * nu_i;

        return value;
    };

    double ***functional_nu_sec_tri_18 = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18[i] = new double *[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18[i][j] = new double [18];
        }
    }

    double ****B_2_tri_3_18 = new double ***[2];
    for (int i = 0; i < 2; i++) {
        B_2_tri_3_18[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            B_2_tri_3_18[i][j] = new double *[3];
            for (int k = 0; k < 3; k++) {
                B_2_tri_3_18[i][j][k] = new double [18]();
            }
        }
    }

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int t = 0; t < 3; t++) {
                for (int k = 0; k < 18*num_modes; k++) {
                    #pragma omp atomic
                    B_2_tri_3_18[0][j][t][k/num_modes] += B_3_tri_18N[t][j][k] *  fourier_basis_test(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
                    #pragma omp atomic
                    B_2_tri_3_18[1][j][t][k/num_modes] += B_3_tri_18N[t][j][k] * fourier_basis_test2(k%num_modes, nfp_ft, (2*pi/num_sections) * i);
                }
            }
        }

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            
            for (int k = 0; k < 18; k++) {
                functional_nu_sec_tri_18[i][j][k] = numerical_integration_divergence_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, G, i, j, num_modes, nfp_ft, num_sections, B_2_tri_3_18, k, \
                                                                                             GLs, n);
            } 
        }

        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int s = 0; s < 2; s++) {
                for (int t = 0; t < 3; t++) {
                    for (int k = 0; k < 18; k++) {
                        B_2_tri_3_18[s][j][t][k] = 0.0;
                    }
                }
            }
        }

    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                delete [] B_2_tri_3_18[i][j][k] ;
            }
            delete [] B_2_tri_3_18[i][j];
        }
        delete [] B_2_tri_3_18[i];
    }
    delete [] B_2_tri_3_18;


    return functional_nu_sec_tri_18;
}

//
double *numerical_integration_poisson_over_triangle(double A, double B, double C, std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int, int, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri, int column, int num_variables, \
double **GLs, int n) {
    // order = 5

    int num = num_variables;

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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, position_sec, position_tri, column, num_variables);
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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, position_sec, position_tri, column, num_variables);
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

double ****integral_poisson_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double *(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double **triangles_uv, int num_variables) {

    printf("calculating integral of (Poisson)*nu_i...\n");
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double*(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int, int, int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri, int column, int num_variables) {

        double nu_i = 0;

        double *value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, position_sec, position_tri);
        for (int i = 0; i < num_variables; i++) {
            value[i] *= nu_i;
        }

        return value;
    };

    double ****functional_nu_sec_tri_18_m = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_m[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_m[i][j] = new double *[18];
        }
    }

    for (int i = 0; i < num_sections; i++) {
        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_m[i][j][k] = numerical_integration_poisson_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, i, j, k, num_variables, \
                                                                                             GLs, n);
            } 
        }
    }

    return functional_nu_sec_tri_18_m;
}

//
double *numerical_integration2_over_triangle(double A, double B, double C, std::function<double *(double, double, double, double, double, double**, double**, double**, int, int, double, double, int, int)> func, \
double x0, double y0, double zeta, double **fc, double **fc_zeta, double **fc_zetazeta, int num_variables, int nfp, double sine, double cosine, int position, int column, \
double **GLs, int n) {
    // order = 5

    int num = num_variables;

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
                ptr = func(u, v, x0, y0, zeta, fc, fc_zeta, fc_zetazeta, num, nfp, sine, cosine, position, column);
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
                ptr = func(u, v, x0, y0, zeta, fc, fc_zeta, fc_zetazeta, num, nfp, sine, cosine, position, column);
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

double ****integral_functional2_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double *(double, double, double, double, double, double **, double **, double **, int, int, double, double)> functional, \
int nfp, double **fc, double **fc_zeta, double **fc_zetazeta, int num_variables, double **triangles_uv) {

    printf("calculating integral of (R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta)*nu_i...\n");
    int num_triangles_in_plane = (int) triangles_uv[0][6];
    int num = num_variables;
    
    std::function<double *(double, double, double, double, double, double**, double **, double **, int, int, double, double, int, int)> functional_nu = \
    [&](double u, double v, double x0, double y0, double zeta, double **fc, double **fc_zeta, double **fc_zetazeta, int num_variables, int nfp, double sine, double cosine, int position, int column) {

        double nu_i = 0;
        double *value;
        for (int j = 0; j < 20; j++) {
            nu_i += (G[position][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, x0, y0, zeta, fc, fc_zeta, fc_zetazeta, num_variables, nfp, sine, cosine);
        for (int i = 0; i < num; i++) {
            value[i] *= nu_i;
        }

        return value;
    }; 

    double ****functional_nu_sec_tri_18_num = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_num[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_num[i][j] = new double *[18];            
        }         
    }

    for (int i = 0; i < num_sections; i++) {
        #pragma omp parallel for num_threads (NT)
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {
                functional_nu_sec_tri_18_num[i][j][k] = numerical_integration2_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                triangles_uv[j][4], triangles_uv[j][5], (double) i * (2*pi / num_sections), fc, fc_zeta, fc_zetazeta, num, nfp, sin(triangles_uv[j][3]), cos(triangles_uv[j][3]), j, k, \
                GLs, n);
            } 
        }
    }

    return functional_nu_sec_tri_18_num;
}

//
double numerical_integration_jacobian_over_triangle(double A, double B, double C, std::function<double(double, double, double ***, double ***, double ***, int, int, int)> func, \
double ***R, double ***Z, double ***G, int position_sec, int position_tri, int column, \
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
                h[k] += func(u, v, R, Z, G, position_sec, position_tri, column);
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
                h[k] += func(u, v, R, Z, G, position_sec, position_tri, column);
            }
        }

        val += c[k] * h[k];  

    }

    val /= (4*n*n);
    val *= (A+B)*C/4;

    return val;
}

double ***integral_jacobian_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double(double u, double v, double ***R, double ***Z, double ***G, int position_sec, int position_tri)> functional, \
double ***R, double ***Z, double **triangles_uv) {

    printf("calculating integral of jacobian*nu_i...\n");
    int num_triangles_in_plane = (int) triangles_uv[0][6];
    
    std::function<double(double , double , double ***, double ***, double ***, int, int, int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***G, int position_sec, int position_tri, int column) {

        double nu_i = 0, value;
        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }
        value = functional(u, v, R, Z, G, position_sec, position_tri) * nu_i;

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
        #pragma omp parallel for num_threads (NT)
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {
                // printf("i = %d, j = %d, k = %d\n", i, j, k);

                functional_nu_sec_tri_18[i][j][k] = numerical_integration_jacobian_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                R, Z, G, i, j, k, \
                GLs, n);
            } 
        }
    }

    return functional_nu_sec_tri_18;
}

//
double ***numerical_integration_beltrami_over_triangle(double A, double B, double C, std::function<double***(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int, int, double, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri, int num_scalar, double mu, \
int column, \
double **GLs, int n) {
    // order = 5

    int num_derivates = 11;

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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, position_sec, position_tri, num_scalar, mu, column);
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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, position_sec, position_tri, num_scalar, mu, column);
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

double ******integral_beltrami_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double ***(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int, int, double)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double **triangles_uv, int num_scalar, double mu) {

    printf("calculating integral of (Beltrami)*nu_i...\n");
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double***(double, double, double ***, double ***, double ***, double ***, double ***, double ***, double ***, int, int, int, double, int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri, int num_scalar, double mu, \
    int column) { 

        int num_derivaties = 11;

        double nu_i = 0;

        double ***value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, position_sec, position_tri, num_scalar, mu);
        for (int i = 0; i < num_derivaties; i++) {
            for (int j = 0; j < num_scalar; j++) {
                for (int k = 0; k < num_scalar; k++) {
                    value[i][j][k] *= nu_i;
                }
            }
            
        }

        return value;
    };

    double ******functional_nu_sec_tri_18_11_3_3 = new double *****[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_11_3_3[i] = new double ****[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_11_3_3[i][j] = new double ***[18];
        }
    }

    for (int i = 0; i < num_sections; i++) {
        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_11_3_3[i][j][k] = numerical_integration_beltrami_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, Rzetazeta, Zzetazeta, G, i, j, num_scalar, mu, \
                                                                                             k, \
                                                                                             GLs, n);
            } 
        }
    }

    return functional_nu_sec_tri_18_11_3_3;
}

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
std::function<double(double u, double v, int position_tr, double **triangles_uv, vector<double> coefficients)> functional, \
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
double ***numerical_integration_SMIEmf_over_triangle(double A, double B, double C, std::function<double***(double , double , double ***, double ***, double ***, double ***, double ***, double **, double ***, int , int , int , int , int , int, int)> func, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***mf_3_tri_18, double **pres_tri_18, double ***G, int num_triangles, int num_modes, int nfp_ft, int position_sec, int position_tri, int num_scalar, \
int column, \
double **GLs, int n) {
    // order = 5

    int num_derivates = 5;

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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, mf_3_tri_18, pres_tri_18, G, num_triangles, num_modes, nfp_ft , position_sec, position_tri, num_scalar, column);
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
                ptr = func(u, v, R, Z, Rzeta, Zzeta, mf_3_tri_18, pres_tri_18, G, num_triangles, num_modes, nfp_ft , position_sec, position_tri, num_scalar, column);
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

double ******integral_SMIEmf_nu_dudv(double ***G, int num_sections, double **GLs, int n, \
std::function<double ***(double , double , double ***, double ***, double ***, double ***, double ***, double **, double ***, int , int , int , int , int , int)> functional, \
double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***mf_3_tri_18N, double **pres_tri_18N, double **triangles_uv, int num_modes, int nfp_ft, int num_scalar) {

    printf("calculating integral of (SMIEmf)*nu_i...\n");
    fflush(stdout);
    int num_triangles_in_plane = (int) triangles_uv[0][6];

    std::function<double***(double , double , double ***, double ***, double ***, double ***, double ***, double **, double ***, int , int , int , int , int , int, int)> functional_nu = \
    [&](double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***mf_3_tri_18, double **pres_tri_18, double ***G, int num_triangles, int num_modes, int nfp_ft, int position_sec, int position_tri, int num_scalar, \
    int column) {

        int num_derivaties = 5;

        double nu_i = 0;

        double ***value;

        for (int j = 0; j < 20; j++) {
            nu_i += (G[position_tri][j][column] * power_nonnegative(u, m_order[j]) * power_nonnegative(v, n_order[j]));
        }

        value = functional(u, v, R, Z, Rzeta, Zzeta, mf_3_tri_18, pres_tri_18, G, num_triangles, num_modes, nfp_ft, position_sec, position_tri, num_scalar);
        for (int i = 0; i < num_derivaties; i++) {
            for (int j = 0; j < num_scalar; j++) {
                for (int k = 0; k < num_scalar; k++) {
                    value[i][j][k] *= nu_i;
                }
            }
            
        }

        return value;
    };

    double ******functional_nu_sec_tri_18_11_3_3 = new double *****[num_sections];
    for (int i = 0; i < num_sections; i++) {
        functional_nu_sec_tri_18_11_3_3[i] = new double ****[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            functional_nu_sec_tri_18_11_3_3[i][j] = new double ***[18];
        }
    }

    double dzeta = 2 * pi / num_sections;
    double *zeta = new double [num_sections]();
    for (int i = 0; i < num_sections; i++) {
        zeta[i] = dzeta * i;
    }

    double ****mf_sec_3_tri_18 = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        mf_sec_3_tri_18[i] = new double **[3];
        for (int j = 0; j < 3; j++) {
            mf_sec_3_tri_18[i][j] = new double *[num_triangles_in_plane];
            for (int k = 0; k < num_triangles_in_plane; k++) {
                mf_sec_3_tri_18[i][j][k] = new double [18]();
            }
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < num_triangles_in_plane; k++) {
                for (int l = 0; l < 18*num_modes; l++) {
                    #pragma omp atomic
                    mf_sec_3_tri_18[i][j][k][l/num_modes] += mf_3_tri_18N[j][k][l] * fourier_basis_test(l%num_modes, nfp_ft, zeta[i]);
                }
            }
        }
    }

    double ***p_sec_tri_18 = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        p_sec_tri_18[i] = new double *[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            p_sec_tri_18[i][j] = new double [18]();
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18*num_modes; k++) {
                #pragma omp atomic
                p_sec_tri_18[i][j][k/num_modes] += pres_tri_18N[j][k] * fourier_basis_test(k%num_modes, nfp_ft, zeta[i]);
            }
        }
    } 

    for (int i = 0; i < num_sections; i++) {

        #pragma omp parallel for
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 18; k++) {

                functional_nu_sec_tri_18_11_3_3[i][j][k] = numerical_integration_SMIEmf_over_triangle(triangles_uv[j][0], triangles_uv[j][1], triangles_uv[j][2], functional_nu, \
                                                                                             R, Z, Rzeta, Zzeta, mf_sec_3_tri_18[i], p_sec_tri_18[i], G, num_triangles_in_plane, num_modes, nfp_ft, i, j, num_scalar, \
                                                                                             k, \
                                                                                             GLs, n);
            } 
        }
    }

    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < num_triangles_in_plane; k++) {
                delete [] mf_sec_3_tri_18[i][j][k];
            }
            delete [] mf_sec_3_tri_18[i][j];
        }
        delete [] mf_sec_3_tri_18[i];
    }
    delete [] mf_sec_3_tri_18;

    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            delete [] p_sec_tri_18[i][j];
        }
        delete [] p_sec_tri_18[i];
    }
    delete [] p_sec_tri_18;

    delete [] zeta;

    return functional_nu_sec_tri_18_11_3_3;
}

