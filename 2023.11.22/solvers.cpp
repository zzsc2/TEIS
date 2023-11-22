#include "solvers.hpp"
#include "superstructure.hpp"

double *FuncPoisson(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri) {
    // extern int m_order[20];
    // extern int n_order[20];

    int dim = 9;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H, alpha[5], beta[5], gamma[5];
    double A_31, A_32, A_41, A_42, A_51, A_52;
    double uv[20], uv_u[20], uv_v[20], uv_uu[20], uv_uv[20], uv_vv[20];
    double R_ = 0, R_u = 0, R_v = 0, R_uu = 0, R_uv = 0, R_vv = 0, R_zetau = 0, R_zetav = 0, R_zetazeta = 0, R_zeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_uu = 0, Z_uv = 0, Z_vv = 0, Z_zetau = 0, Z_zetav = 0, Z_zetazeta = 0, Z_zeta = 0;

    for (int i = 0; i < 20; i++) {
        uv[i] =                                   power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] *                   power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] *                   power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);
        uv_uu[i] = m_order[i] * (m_order[i]-1) *  power_nonnegative(u, m_order[i]-2) * power_nonnegative(v, n_order[i]);
        uv_uv[i] = m_order[i] *  n_order[i] *     power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]-1);
        uv_vv[i] = n_order[i] * (n_order[i]-1) *  power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-2);
        for (int j = 0; j < 18; j++) {
            R_ +=         G[position_tri][i][j] * uv[i] *    R[position_sec][position_tri][j];
            R_u +=        G[position_tri][i][j] * uv_u[i] *  R[position_sec][position_tri][j];
            R_v +=        G[position_tri][i][j] * uv_v[i] *  R[position_sec][position_tri][j];
            R_uu +=       G[position_tri][i][j] * uv_uu[i] * R[position_sec][position_tri][j];
            R_uv +=       G[position_tri][i][j] * uv_uv[i] * R[position_sec][position_tri][j];
            R_vv +=       G[position_tri][i][j] * uv_vv[i] * R[position_sec][position_tri][j];
            
            Z_ +=         G[position_tri][i][j] * uv[i] *    Z[position_sec][position_tri][j];
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

    //R  Z  RR RZ ZZ φR φZ φφ φ
    //0  1  2  3  4  5  6  7  8
    double jacobian =  D * R_;
    double *var1 = new double [10]();
    for (int i = 1; i < 10; i++) {
        var1[i] = (p[0][i-1] / R_ + p[2][i-1] + p[4][i-1] + p[7][i-1] / (R_*R_)) * jacobian;
    }

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}

double ***FuncBeltrami(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, int position_sec, int position_tri, int num_scalar, double mu) {
    // extern int m_order[20];
    // extern int n_order[20];

    int dim = 9;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H;
    double alpha[5], beta[5], gamma[5];
    double A_31, A_32, A_41, A_42, A_51, A_52;
    double uv[20], uv_u[20], uv_v[20];
    double uv_uu[20], uv_uv[20], uv_vv[20];
    double R_ = 0, R_u = 0, R_v = 0, R_zeta = 0;
    double R_uu = 0, R_uv = 0, R_vv = 0, R_zetau = 0, R_zetav = 0, R_zetazeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_zeta = 0;
    double Z_uu = 0, Z_uv = 0, Z_vv = 0, Z_zetau = 0, Z_zetav = 0, Z_zetazeta = 0;

    for (int i = 0; i < 20; i++) {
        uv[i] =                                   power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] *                   power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] *                   power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);
        uv_uu[i] = m_order[i] * (m_order[i]-1) *  power_nonnegative(u, m_order[i]-2) * power_nonnegative(v, n_order[i]);
        uv_uv[i] = m_order[i] *  n_order[i] *     power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]-1);
        uv_vv[i] = n_order[i] * (n_order[i]-1) *  power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-2);
        for (int j = 0; j < 18; j++) {
            R_ +=         G[position_tri][i][j] * uv[i] *    R[position_sec][position_tri][j];
            R_u +=        G[position_tri][i][j] * uv_u[i] *  R[position_sec][position_tri][j];
            R_v +=        G[position_tri][i][j] * uv_v[i] *  R[position_sec][position_tri][j];
            R_uu +=       G[position_tri][i][j] * uv_uu[i] * R[position_sec][position_tri][j];
            R_uv +=       G[position_tri][i][j] * uv_uv[i] * R[position_sec][position_tri][j];
            R_vv +=       G[position_tri][i][j] * uv_vv[i] * R[position_sec][position_tri][j];
            
            Z_ +=         G[position_tri][i][j] * uv[i] *    Z[position_sec][position_tri][j];
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
    // - R_vv*Z_u^2*Z_v + R_v*Z_vv*Z_u^2 + 2*R_uv*Z_u*Z_v^2 - 2*R_v*Z_uv*Z_u*Z_v - R_uu*Z_v^3 + R_v*Z_uu*Z_v^2
    // R_vv*Z_u^3 - 2*R_uv*Z_u^2*Z_v - R_u*Z_vv*Z_u^2 + R_uu*Z_u*Z_v^2 + 2*R_u*Z_uv*Z_u*Z_v - R_u*Z_uu*Z_v^2
    // R_uu*R_v*Z_v^2 - R_u*R_uv*Z_v^2 + R_v^2*Z_u*Z_uv - R_v^2*Z_uu*Z_v - R_u*R_v*Z_u*Z_vv + R_u*R_v*Z_uv*Z_v + R_u*R_vv*Z_u*Z_v - R_uv*R_v*Z_u*Z_v
    // R_uv*R_v*Z_u^2 - R_u*R_vv*Z_u^2 + R_u^2*Z_u*Z_vv - R_u^2*Z_uv*Z_v + R_u*R_uv*Z_u*Z_v - R_u*R_v*Z_u*Z_uv + R_u*R_v*Z_uu*Z_v - R_uu*R_v*Z_u*Z_v
    // Z_vv*R_u^2*R_v - R_vv*Z_v*R_u^2 - 2*Z_uv*R_u*R_v^2 + 2*R_uv*Z_v*R_u*R_v + Z_uu*R_v^3 - R_uu*Z_v*R_v^2
    // - Z_vv*R_u^3 + 2*Z_uv*R_u^2*R_v + R_vv*Z_u*R_u^2 - Z_uu*R_u*R_v^2 - 2*R_uv*Z_u*R_u*R_v + R_uu*Z_u*R_v^2

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

    //R  Z  RR RZ ZZ φR φZ φφ φ
    //0  1  2  3  4  5  6  7  8

    // num_scalar = 3 for beltrami field
    double ***var1 = new double **[11];
    for (int i = 0; i < 11; i++) {
        var1[i] = new double *[num_scalar];
        for (int j = 0; j < num_scalar; j++) {
            var1[i][j] = new double [num_scalar]();
        }
    }


    var1[0][0][0] = (- 1.0 / (R_*R_)) * D * R_;
    var1[0][2][2] = (- 1.0 / (R_*R_)) * D * R_;
    for (int i = 1; i < 10; i++) {
       var1[i][0][0] = (p[0][i-1] / R_ + p[2][i-1] + p[4][i-1] + p[7][i-1] / (R_*R_)) * D * R_;
       var1[i][0][2] = (- 2.0 / (R_*R_) * p[8][i-1]) * D * R_;

       var1[i][1][1] = (p[0][i-1] / R_ + p[2][i-1] + p[4][i-1] + p[7][i-1] / (R_*R_)) * D * R_;

       var1[i][2][0] = (+ 2.0 / (R_*R_) * p[8][i-1]) * D * R_;
       var1[i][2][2] = (p[0][i-1] / R_ + p[2][i-1] + p[4][i-1] + p[7][i-1] / (R_*R_)) * D * R_;
    }
    var1[0][1][2] += (mu / R_) * D * R_;
    for (int i = 1; i < 10; i++) {
        var1[i][0][1] += (p[8][i-1] / R_ * mu) * D * R_;
        var1[i][0][2] += (- p[1][i-1] * mu) * D * R_;

        var1[i][1][0] += (- p[8][i-1] / R_ * mu) * D * R_;
        var1[i][1][2] += (p[0][i-1] * mu) * D * R_;

        var1[i][2][0] += (p[1][i-1] * mu) * D * R_;
        var1[i][2][1] += (- p[0][i-1] * mu) * D * R_;
    }





    // var1[0][1][2] += (1.0 / R_) * D * R_;
    // for (int i = 1; i < 10; i++) {
    //     var1[i][0][1] += (  p[8][i-1] / R_) * D * R_;
    //     var1[i][0][2] += (- p[1][i-1]) * D * R_;

    //     var1[i][1][0] += (- p[8][i-1] / R_) * D * R_;
    //     var1[i][1][2] += (  p[0][i-1]) * D * R_;

    //     var1[i][2][0] += (  p[1][i-1]) * D * R_;
    //     var1[i][2][1] += (- p[0][i-1]) * D * R_;
    // }
    // var1[0][0][0] += (- mu) * D * R_;
    // var1[0][1][1] += (- mu) * D * R_;
    // var1[0][2][2] += (- mu) * D * R_;    







    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}

double FuncDivergence(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ****B_2_tri_3_18) {

    int dim = 3;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H;
    double uv[20], uv_u[20], uv_v[20];
    double R_ = 0, R_u = 0, R_v = 0, R_zeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_zeta = 0;

    // double B_3_18[3][18] = {0};
    // double B_zeta_3_18[3][18] = {0};
    double BR = 0, BR_u = 0, BR_v = 0, BR_zeta = 0;
    double BZ_u = 0, BZ_v = 0, BZ_zeta = 0;
    double Bzeta_u = 0, Bzeta_v = 0, Bzeta_zeta = 0;

    // for (int k = 0; k < 3; k++) {
    //     for (int j = 0; j < 18*num_modes; j++) {
    //         B_3_18[k][j/num_modes] +=      B_3_tri_18N[k][position_tri][j] * fourier_basis_test(j%num_modes, nfp_ft, (2*pi/num_sections) * position_sec);
    //         B_zeta_3_18[k][j/num_modes] += B_3_tri_18N[k][position_tri][j] * fourier_basis_test2(j%num_modes, nfp_ft, (2*pi/num_sections) * position_sec);
    //     }
    // }

    for (int i = 0; i < 20; i++) {
        uv[i] =                 power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] * power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] * power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);

        for (int j = 0; j < 18; j++) {
            R_ +=         G[position_tri][i][j] * uv[i] *    R[position_sec][position_tri][j];
            R_u +=        G[position_tri][i][j] * uv_u[i] *  R[position_sec][position_tri][j];
            R_v +=        G[position_tri][i][j] * uv_v[i] *  R[position_sec][position_tri][j];
            
            Z_ +=         G[position_tri][i][j] * uv[i] *    Z[position_sec][position_tri][j];
            Z_u +=        G[position_tri][i][j] * uv_u[i] *  Z[position_sec][position_tri][j];
            Z_v +=        G[position_tri][i][j] * uv_v[i] *  Z[position_sec][position_tri][j];

            R_zeta +=     G[position_tri][i][j] * uv[i] *    Rzeta[position_sec][position_tri][j];

            Z_zeta +=     G[position_tri][i][j] * uv[i] *    Zzeta[position_sec][position_tri][j];

            // BR +=         G[position_tri][i][j] * uv[i] *    B_3_18[0][j];
            // BR_u +=       G[position_tri][i][j] * uv_u[i] *  B_3_18[0][j];
            // BR_v +=       G[position_tri][i][j] * uv_v[i] *  B_3_18[0][j];
            // BR_zeta +=    G[position_tri][i][j] * uv[i] *    B_zeta_3_18[0][j];

            // BZ_u +=       G[position_tri][i][j] * uv_u[i] *  B_3_18[1][j];
            // BZ_v +=       G[position_tri][i][j] * uv_v[i] *  B_3_18[1][j];
            // BZ_zeta +=    G[position_tri][i][j] * uv[i] *    B_zeta_3_18[1][j];

            // Bzeta_u +=    G[position_tri][i][j] * uv_u[i] *  B_3_18[2][j];
            // Bzeta_v +=    G[position_tri][i][j] * uv_v[i] *  B_3_18[2][j];
            // Bzeta_zeta += G[position_tri][i][j] * uv[i] *    B_zeta_3_18[2][j];

            BR +=         G[position_tri][i][j] * uv[i] *    B_2_tri_3_18[0][position_tri][0][j];

            BR_u +=       G[position_tri][i][j] * uv_u[i] *  B_2_tri_3_18[0][position_tri][0][j];
            BR_v +=       G[position_tri][i][j] * uv_v[i] *  B_2_tri_3_18[0][position_tri][0][j];
            BR_zeta +=    G[position_tri][i][j] * uv[i] *    B_2_tri_3_18[1][position_tri][0][j];

            BZ_u +=       G[position_tri][i][j] * uv_u[i] *  B_2_tri_3_18[0][position_tri][1][j];
            BZ_v +=       G[position_tri][i][j] * uv_v[i] *  B_2_tri_3_18[0][position_tri][1][j];
            BZ_zeta +=    G[position_tri][i][j] * uv[i] *    B_2_tri_3_18[1][position_tri][1][j];

            Bzeta_u +=    G[position_tri][i][j] * uv_u[i] *  B_2_tri_3_18[0][position_tri][2][j];
            Bzeta_v +=    G[position_tri][i][j] * uv_v[i] *  B_2_tri_3_18[0][position_tri][2][j];
            Bzeta_zeta += G[position_tri][i][j] * uv[i] *    B_2_tri_3_18[1][position_tri][2][j];
            
        }
    }

    D = R_u * Z_v - R_v * Z_u;
    E = R_v * Z_zeta - R_zeta * Z_v;
    H = R_zeta * Z_u - R_u * Z_zeta;

    p[0][0] = Z_v / D;
    p[0][1] = - Z_u / D;

    p[1][0] = - R_v / D;
    p[1][1] = R_u / D;


    p[2][0] = E / D;
    p[2][1] = H / D;
    p[2][2] = 1;

    //R  Z  φ
    //0  1  2

    double var1 = 0.0;

    var1 += 1.0 / R_ * BR;

    var1 += p[0][0] * BR_u;
    var1 += p[0][1] * BR_v;
    var1 += p[0][2] * BR_zeta;

    var1 += p[1][0] * BZ_u;
    var1 += p[1][1] * BZ_v;
    var1 += p[1][2] * BZ_zeta;

    var1 += p[2][0] / R_ * Bzeta_u;
    var1 += p[2][1] / R_ * Bzeta_v;
    var1 += p[2][2] / R_ * Bzeta_zeta;

    var1 *= (D * R_);

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}

double *FuncProjection(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ***B_tri_3_18, double ***poisoon_2_tri_18) {

    int dim = 3;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H;
    double uv[20], uv_u[20], uv_v[20];
    double R_ = 0, R_u = 0, R_v = 0, R_zeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_zeta = 0;

    // double B_3_18[3][18] = {0};
    // double poisson_18[18] = {0};
    // double poisson_zeta_18[18] = {0};
    double BR = 0, BZ = 0, Bzeta = 0;
    double poisson_u = 0, poisson_v = 0, poisson_zeta = 0;

    // for (int k = 0; k < 3; k++) {
    //     for (int j = 0; j < 18*num_modes; j++) {
    //         B_3_18[k][j/num_modes] += B_3_tri_18N[k][position_tri][j] * fourier_basis_test(j%num_modes, nfp_ft, (2*pi/num_sections) * position_sec);
    //     }
    // }
    // for (int j = 0; j < 18*num_modes; j++) {
    //     poisson_18[j/num_modes] +=      poisoon_1_tri_18N[0][position_tri][j] * fourier_basis_test(j%num_modes, nfp_ft, (2*pi/num_sections) * position_sec);
    //     poisson_zeta_18[j/num_modes] += poisoon_1_tri_18N[0][position_tri][j] * fourier_basis_test2(j%num_modes, nfp_ft, (2*pi/num_sections) * position_sec);
    // }

    for (int i = 0; i < 20; i++) {
        uv[i] =                 power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] * power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] * power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);

        for (int j = 0; j < 18; j++) {
            R_ +=         G[position_tri][i][j] * uv[i] *    R[position_sec][position_tri][j];
            R_u +=        G[position_tri][i][j] * uv_u[i] *  R[position_sec][position_tri][j];
            R_v +=        G[position_tri][i][j] * uv_v[i] *  R[position_sec][position_tri][j];
            
            Z_ +=         G[position_tri][i][j] * uv[i] *    Z[position_sec][position_tri][j];
            Z_u +=        G[position_tri][i][j] * uv_u[i] *  Z[position_sec][position_tri][j];
            Z_v +=        G[position_tri][i][j] * uv_v[i] *  Z[position_sec][position_tri][j];

            R_zeta +=     G[position_tri][i][j] * uv[i] *    Rzeta[position_sec][position_tri][j];

            Z_zeta +=     G[position_tri][i][j] * uv[i] *    Zzeta[position_sec][position_tri][j];

            BR +=         G[position_tri][i][j] * uv[i] *    B_tri_3_18[position_tri][0][j];
            BZ +=         G[position_tri][i][j] * uv[i] *    B_tri_3_18[position_tri][1][j];
            Bzeta +=      G[position_tri][i][j] * uv[i] *    B_tri_3_18[position_tri][2][j];

            poisson_u +=  G[position_tri][i][j] * uv_u[i] *  poisoon_2_tri_18[0][position_tri][j];
            poisson_v +=  G[position_tri][i][j] * uv_v[i] *  poisoon_2_tri_18[0][position_tri][j];
          poisson_zeta += G[position_tri][i][j] * uv[i] *    poisoon_2_tri_18[1][position_tri][j];
        //     poisson_u +=  G[position_tri][i][j] * uv_u[i] *  poisson_18[j];
        //     poisson_v +=  G[position_tri][i][j] * uv_v[i] *  poisson_18[j];
        //   poisson_zeta += G[position_tri][i][j] * uv[i] *    poisson_zeta_18[j];

        }
    }

    D = R_u * Z_v - R_v * Z_u;
    E = R_v * Z_zeta - R_zeta * Z_v;
    H = R_zeta * Z_u - R_u * Z_zeta;

    p[0][0] = Z_v / D;
    p[0][1] = - Z_u / D;

    p[1][0] = - R_v / D;
    p[1][1] = R_u / D;


    p[2][0] = E / D;
    p[2][1] = H / D;
    p[2][2] = 1;

    //R  Z  φ
    //0  1  2

    double *var1 = new double [3]();

    var1[0] = (BR -    (poisson_u * p[0][0] + poisson_v * p[0][1]));
    var1[1] = (BZ -    (poisson_u * p[1][0] + poisson_v * p[1][1]));
    var1[2] = (Bzeta - (poisson_u * p[2][0] + poisson_v * p[2][1] + poisson_zeta) / R_);  

    // var1[0] = (BR -    (0));
    // var1[1] = (BZ -    (0));
    // var1[2] = (Bzeta - (0) / R_); 

    // var1[0] = (0 -    (1 * poisson_u + 1 * poisson_u + 1 * poisson_u));
    // var1[1] = (0 -    (1 * poisson_v + 1 * poisson_v + 1 * poisson_v));
    // var1[2] = (0 -    (1 * poisson_zeta + 1 * poisson_zeta + 1 * poisson_zeta) / R_);

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}

double ***FuncProjection2(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, int position_sec, int position_tri, int num_modes, int nfp_ft, int num_sections, double ***B_tri_3_18, double ***poisoon_2_tri_18) {

    int dim = 3;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H;
    double uv[20], uv_u[20], uv_v[20];
    double R_ = 0, R_u = 0, R_v = 0, R_zeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_zeta = 0;

    double BR = 0, BZ = 0, Bzeta = 0;
    double poisson_u = 0, poisson_v = 0, poisson_zeta = 0;

    for (int i = 0; i < 20; i++) {
        uv[i] =                 power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] * power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] * power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);

        for (int j = 0; j < 18; j++) {
            R_ +=         G[position_tri][i][j] * uv[i] *    R[position_sec][position_tri][j];
            R_u +=        G[position_tri][i][j] * uv_u[i] *  R[position_sec][position_tri][j];
            R_v +=        G[position_tri][i][j] * uv_v[i] *  R[position_sec][position_tri][j];
            
            Z_ +=         G[position_tri][i][j] * uv[i] *    Z[position_sec][position_tri][j];
            Z_u +=        G[position_tri][i][j] * uv_u[i] *  Z[position_sec][position_tri][j];
            Z_v +=        G[position_tri][i][j] * uv_v[i] *  Z[position_sec][position_tri][j];

            R_zeta +=     G[position_tri][i][j] * uv[i] *    Rzeta[position_sec][position_tri][j];

            Z_zeta +=     G[position_tri][i][j] * uv[i] *    Zzeta[position_sec][position_tri][j];

            BR +=         G[position_tri][i][j] * uv[i] *    B_tri_3_18[position_tri][0][j];
            BZ +=         G[position_tri][i][j] * uv[i] *    B_tri_3_18[position_tri][1][j];
            Bzeta +=      G[position_tri][i][j] * uv[i] *    B_tri_3_18[position_tri][2][j];

            poisson_u +=  G[position_tri][i][j] * uv_u[i] *  poisoon_2_tri_18[0][position_tri][j];
            poisson_v +=  G[position_tri][i][j] * uv_v[i] *  poisoon_2_tri_18[0][position_tri][j];
          poisson_zeta += G[position_tri][i][j] * uv[i] *    poisoon_2_tri_18[1][position_tri][j];

        }
    }

    D = R_u * Z_v - R_v * Z_u;
    E = R_v * Z_zeta - R_zeta * Z_v;
    H = R_zeta * Z_u - R_u * Z_zeta;

    p[0][0] = Z_v / D;
    p[0][1] = - Z_u / D;

    p[1][0] = - R_v / D;
    p[1][1] = R_u / D;


    p[2][0] = E / D;
    p[2][1] = H / D;
    p[2][2] = 1;

    //R  Z  φ
    //0  1  2

    double ***var1 = new double **[2];
    for (int i = 0; i < 2; i++) {
        var1[i] = new double *[3];
        for (int j = 0; j < 3; j++) {
            var1[i][j] = new double [3]();
        }
    }


    var1[0][0][0] = 1.0 * D * R_;
    var1[0][1][1] = 1.0 * D * R_;
    var1[0][2][2] = 1.0 * D * R_;  

    var1[1][0][0] = (BR -    (poisson_u * p[0][0] + poisson_v * p[0][1])) * D * R_;
    var1[1][1][0] = (BZ -    (poisson_u * p[1][0] + poisson_v * p[1][1])) * D * R_;
    var1[1][2][0] = (Bzeta - (poisson_u * p[2][0] + poisson_v * p[2][1] + poisson_zeta) / R_) * D * R_;  

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}

double *FuncDivCoeff(int position_p, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, int position_sec, int position_tri, int num_sections) {

    int dim = 3;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H;
    double R_, R_x, R_y, R_zeta;
    double Z_, Z_x, Z_y, Z_zeta;

    R_ =         R[position_sec][position_tri][6*position_p+0];
    R_x =        R[position_sec][position_tri][6*position_p+1];
    R_y =        R[position_sec][position_tri][6*position_p+2];

    Z_ =         Z[position_sec][position_tri][6*position_p+0];
    Z_x =        Z[position_sec][position_tri][6*position_p+1];
    Z_y =        Z[position_sec][position_tri][6*position_p+2];

    R_zeta = Rzeta[position_sec][position_tri][6*position_p+0];

    Z_zeta = Zzeta[position_sec][position_tri][6*position_p+0];

    E = R_y * Z_zeta - R_zeta * Z_y;
    H = R_zeta * Z_x - R_x * Z_zeta;
    D = R_x * Z_y - R_y * Z_x;

    p[0][0] =   Z_y / D;
    p[0][1] = - Z_x / D;

    p[1][0] = - R_y / D;
    p[1][1] =   R_x / D;


    p[2][0] = E / D;
    p[2][1] = H / D;
    p[2][2] = 1.0;

    //R  Z  φ
    //0  1  2

    double *var1 = new double [9]();

    var1[0] = 1.0 / R_;
    var1[1] = 0.0;
    var1[2] = p[2][2] / R_;

    var1[3] = p[0][0];
    var1[4] = p[1][0];
    var1[5] = p[2][0] / R_;

    var1[6] = p[0][1];
    var1[7] = p[1][1];
    var1[8] = p[2][1] / R_;

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}

double ****divergence_free_constraints(double **triangles_uv, int num_sections, int num_modes, int nfp_ft, double ***R, double ***Z, double ***Rzeta, double ***Zzeta) {
    printf("constructing divergence-free constraints...\n");

    int num_triangles_in_plane = (int) triangles_uv[0][6];
    int num_variables = 9;
    int num_modesD = 2 * num_modes - 1;

    double ****div_sec_tri_3_var = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        div_sec_tri_3_var[i] = new double **[num_triangles_in_plane];
        for (int j = 0; j < num_triangles_in_plane; j++) {
            div_sec_tri_3_var[i][j] = new double *[3];
        }
    }

    double ****div_tri_3_var_modes = new double ***[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        div_tri_3_var_modes[i] = new double **[3];
        for (int j = 0; j < 3; j++) {
            div_tri_3_var_modes[i][j] = new double *[num_variables];
            for (int k = 0; k < num_variables; k++) {
                div_tri_3_var_modes[i][j][k] = new double [num_modesD]();
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                div_sec_tri_3_var[i][j][k] = FuncDivCoeff(k, R, Z, Rzeta, Zzeta, i, j, num_sections);
            }
        }
    }
    
    double dzeta = 2 * pi / num_sections;
    double *zeta = new double [num_sections]();
    for (int i = 0; i < num_sections; i++) {
        zeta[i] = dzeta * i;
    }
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < num_triangles_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < num_variables; k++) {
                for (int s = 0; s < num_modesD; s++) {
                    for (int t = 0; t < num_sections; t++) {
                        #pragma omp atomic
                        div_tri_3_var_modes[i][j][k][s] += div_sec_tri_3_var[t][i][j][k] * fourier_basis(s, nfp_ft, zeta[t]) * dzeta;
                    }
                }
            }
        }
    }

    for (int i = 0; i < num_sections; i++) {
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for (int k = 0; k < 3; k++) {
                delete [] div_sec_tri_3_var[i][j][k];
            }
            delete [] div_sec_tri_3_var[i][j];
        }
        delete [] div_sec_tri_3_var[i];
    }
    delete [] div_sec_tri_3_var;

    delete [] zeta;

    // for (int i = 0; i < num_triangles_in_plane; i++) {
    //     for (int j = 0; j < num_variables; j++) {
    //         for (int k = 0; k < num_modesD; k++) {
    //             printf("%0.9f ", div_tri_3_var_modes[i][0][j][k]);
    //         }
    //         printf("   ");
    //     }
    //     printf("\n");
    // }
    // fflush(stdout);
    // exit(0);

    return div_tri_3_var_modes;
}





double **boundary_DoFs_beltrami(double ***triangles_RZ, int **points_sequence, double ***coils_p, int num_sections, int num_modes, int nfp_ft) {
    printf("constructing boundary conditions...\n");

    int num_edges = points_sequence[0][9];
    int num_vertices = points_sequence[0][10];
    int num_triangles_in_plane = (int) triangles_RZ[0][0][6];
    int num_fixed_DoFs = num_edges*num_modes;
    double dzeta = 2 * pi / num_sections;
    double zeta;

    double **value_3_fixed = new double *[3];
    double ***fields_3_sec_edges = new double **[3];
    for (int i = 0; i < 3; i++) {
        value_3_fixed[i] = new double [num_fixed_DoFs]();
        fields_3_sec_edges[i] = new double *[num_sections];
        for (int j = 0; j < num_sections; j++) {
            fields_3_sec_edges[i][j] = new double [num_edges]();
        }
    }

    double *temp1;
    for (int i = 0; i < num_sections; i++) {
        zeta = i * dzeta;
        for (int j = 0; j < num_triangles_in_plane; j++) {
            for(int k = 0; k < 3; k++) {
                if (points_sequence[j][0+3*k] == 1) {
                    temp1 = biot_savart(triangles_RZ[i][j][0+2*k], triangles_RZ[i][j][1+2*k], zeta, coils_p);
                    fields_3_sec_edges[0][i][points_sequence[j][1+3*k]] = temp1[0];
                    fields_3_sec_edges[1][i][points_sequence[j][1+3*k]] = temp1[1];
                    fields_3_sec_edges[2][i][points_sequence[j][1+3*k]] = temp1[2];
                    delete [] temp1;
                }
            }
        }
    }

    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < num_edges; j++) {
    //         for (int k = 0 ; k < num_sections; k++) {
    //             printf("%0.9f  ", fields_3_sec_edges[i][k][j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n\n\n");
    // }
    // exit(0);

    for (int i = 0; i < num_fixed_DoFs; i++) {
        for (int j = 0; j < num_sections; j++) {
            zeta = dzeta * j;
            value_3_fixed[0][i] += fields_3_sec_edges[0][j][i/num_modes] * fourier_basis(i%num_modes, nfp_ft, zeta) * dzeta;
            value_3_fixed[1][i] += fields_3_sec_edges[1][j][i/num_modes] * fourier_basis(i%num_modes, nfp_ft, zeta) * dzeta;
            value_3_fixed[2][i] += fields_3_sec_edges[2][j][i/num_modes] * fourier_basis(i%num_modes, nfp_ft, zeta) * dzeta;
        }
    }

    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < num_fixed_DoFs; j++) {
    //         printf("%0.9f ", value_3_fixed[i][j]);
    //         if ((j+1) % num_modes == 0) {
    //             printf("\n");
    //         }
    //     }
    //     printf("\n\n\n");
    // }
    // exit(0);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < num_sections; j++) {
            delete [] fields_3_sec_edges[i][j];
        }
        delete [] fields_3_sec_edges[i];
    }
    delete [] fields_3_sec_edges;

    return value_3_fixed;
}

double *****MG_tri_6_18_18_18(double ***G, double **triangles_uv) {
    printf("calculating MG...\n");
    fflush(stdout); 

    int num_triangles_in_plane = (int) triangles_uv[0][6];
    double *****MG = new double ****[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        MG[i] = new double ***[6];
        for (int j = 0; j < 6; j++) {
            MG[i][j] = new double **[18];
            for (int k = 0; k < 18; k++) {
                MG[i][j][k] = new double *[18];
                for (int l = 0; l < 18; l++) {
                    MG[i][j][k][l] = new double [18]();
                }
            }
        }
    }

    double *****array_integral3 = new double ****[20];
    for (int i = 0; i < 20; i++) {
        array_integral3[i] = new double ***[20];
        for (int j = 0; j < 20; j++) {
            array_integral3[i][j] = new double **[20];
            for (int k = 0; k < 20; k++) {
                array_integral3[i][j][k] = new double *[3];
                for (int l = 0; l < 3; l++) {
                    array_integral3[i][j][k][l] = new double [3]();
                }
            }
        }
    }
   
    for (int i = 0; i < num_triangles_in_plane; i++) {

        for (int l = 0; l < 20; l++) {
            for (int s = 0; s < 3; s++) {
                for (int t = 0; t < 3; t++) {
                    if ((m_order[l] >= s) && (n_order[l] >= t)) {
                        #pragma omp parallel for collapse(2)
                        for (int j = 0; j < 20; j++) {
                            for (int k = 0; k < 20; k++) {
                                array_integral3[j][k][l][s][t] = \
                                integral1(m_order[j] + m_order[k] + m_order[l] - s, n_order[j] + n_order[k] + n_order[l] - t, triangles_uv[i][0], triangles_uv[i][1], triangles_uv[i][2]);
                            }
                        }
                    }
                }
            }
        }

        #pragma omp parallel for collapse(3)
        for (int j = 0; j < 18; j++) {
            for (int k = 0; k < 18; k++) {
                for (int l = 0; l < 18; l++) {
                    MG[i][0][j][k][l] = integral3(G, j, k, l, i,  0,  0, array_integral3);
                    MG[i][1][j][k][l] = integral3(G, j, k, l, i, -1,  0, array_integral3);
                    MG[i][2][j][k][l] = integral3(G, j, k, l, i,  0, -1, array_integral3);
                    MG[i][3][j][k][l] = integral3(G, j, k, l, i, -2,  0, array_integral3);
                    MG[i][4][j][k][l] = integral3(G, j, k, l, i, -1, -1, array_integral3);
                    MG[i][5][j][k][l] = integral3(G, j, k, l, i,  0, -2, array_integral3);
                }
            }
        }

    }

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            for (int k = 0; k < 20; k++) {
                for (int l = 0; l < 3; l++) {
                    delete [] array_integral3[i][j][k][l];
                }
                delete [] array_integral3[i][j][k];
            }
            delete [] array_integral3[i][j];
        }
        delete [] array_integral3[i];
    }
    delete [] array_integral3;

    return MG;
}

double ***BG_tri_18_18(double ***G, double **triangles_uv) {
    printf("calculating BG...\n");
    fflush(stdout); 
    // extern double integral2(double a, double b, double c, double ***G, int i, int j, int position);

    int num_triangles_in_plane = (int) triangles_uv[0][6];
    double ***B = new double **[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        B[i] = new double *[18];
        for (int j = 0; j < 18; j++) {
            B[i][j] = new double [18]();
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < num_triangles_in_plane; i++) {
        for (int j = 0; j < 18; j++) {
            for (int k = 0; k < 18; k++) {
                B[i][j][k] = integral2(triangles_uv[i][0], triangles_uv[i][1], triangles_uv[i][2], G, j, k, i);
            }
        }
    }
    
    // for (int i = 0; i < num_triangles_in_plane; i++) {
    //     for (int j = 0; j < 18; j++) {
    //         for (int k = 0; k < 18; k++) {
    //             printf("%0.9f ", B[i][j][k]);
    //         }
    //     }
    //     printf("\n");
    // }

    return B;
}

double ****M_integral4(int nfp, int num_modes, double error) {
    printf("calculating M2_integral4...\n");
    fflush(stdout);
    double ****array_integral4 = new double ***[18*num_modes];
    for (int i = 0; i < 18*num_modes; i++) {
        array_integral4[i] = new double **[18*num_modes];
        for (int j = 0; j < 18*num_modes; j++) {
            array_integral4[i][j] = new double *[18*num_modes];
            for (int k = 0; k < 18*num_modes; k++) {
                array_integral4[i][j][k] = new double [3];
                for (int l = 0; l < 3; l++) {
                    array_integral4[i][j][k][l] = integral4(nfp, num_modes, i, j, k, l, error);
                }
            } 
        }
    }
    return array_integral4;
}

double ***M_tri_18Nn_18Nn(double *****MG_tri_6_18_18_18, double *****object_tri_18N_ope_3_3, double ****array_integral4, int num_triangles_in_plane, int num_modes, int nfp, int num_scalar, double error, int num_derivates) {
    printf("calculating M2...\n");
    fflush(stdout); 

    double ***M = new double **[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        M[i] = new double *[(18*num_modes)*num_scalar];
        for (int k = 0; k < (18*num_modes)*num_scalar; k++) {
            M[i][k] = new double [(18*num_modes)*num_scalar]();
        }
    }

    int *tkj = new int [18*num_modes];
    for (int i = 0; i < 18*num_modes; i++) {
        tkj[i] = i/num_modes;
    }

    if (num_derivates == 10) {
        #pragma omp parallel for
        for (int i = 0; i < num_triangles_in_plane; i++) {
            for (int t = 0; t < 18*num_modes; t++) {
                for (int k = 0; k < 18*num_modes; k++) {
                     for (int j = 0; j < 18*num_modes; j++){
                        for (int p = 0; p < num_scalar; p++) {
                            for (int q = 0; q < num_scalar; q++) {
                                // u  v  uu uv vv ζu ζv ζζ ζ
                                // 1  2  3  4  5  6  7  8  9
                                #pragma omp atomic
                                M[i][k*num_scalar + p][j*num_scalar + q] += \
                                object_tri_18N_ope_3_3[i][t][0][p][q] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][0][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][1][p][q] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][1][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][2][p][q] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][2][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][3][p][q] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][3][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][4][p][q] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][4][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][5][p][q] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][5][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][9][p][q] * array_integral4[t][k][j][1] * MG_tri_6_18_18_18[i][0][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][6][p][q] * array_integral4[t][k][j][1] * MG_tri_6_18_18_18[i][1][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][7][p][q] * array_integral4[t][k][j][1] * MG_tri_6_18_18_18[i][2][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][8][p][q] * array_integral4[t][k][j][2] * MG_tri_6_18_18_18[i][0][tkj[t]][tkj[k]][tkj[j]];
                            }
                        }
                    }
                }
            }
        }
    }
    else if (num_derivates == 4) {
        #pragma omp parallel for
        for (int i = 0; i < num_triangles_in_plane; i++) {
            for (int t = 0; t < 18*num_modes; t++) {
                for (int k = 0; k < 18*num_modes; k++) {
                     for (int j = 0; j < 18*num_modes; j++){
                        for (int p = 0; p < num_scalar; p++) {
                            for (int q = 0; q < num_scalar; q++) {
                                #pragma omp atomic
                                M[i][k*num_scalar + p][j*num_scalar + q] += \
                                object_tri_18N_ope_3_3[i][t][0][p][q] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][0][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][1][p][q] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][1][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][2][p][q] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][2][tkj[t]][tkj[k]][tkj[j]] + \
                                object_tri_18N_ope_3_3[i][t][3][p][q] * array_integral4[t][k][j][1] * MG_tri_6_18_18_18[i][0][tkj[t]][tkj[k]][tkj[j]];
                            }
                        }
                    }
                }
            }
        }
    }
    else if (num_derivates == 1) {
        #pragma omp parallel for
        for (int i = 0; i < num_triangles_in_plane; i++) {
            for (int t = 0; t < 18*num_modes; t++) {
                for (int k = 0; k < 18*num_modes; k++) {
                     for (int j = 0; j < 18*num_modes; j++){
                        for (int p = 0; p < num_scalar; p++) {
                            for (int q = 0; q < num_scalar; q++) {
                                #pragma omp atomic
                                M[i][k*num_scalar + p][j*num_scalar + q] += \
                                object_tri_18N_ope_3_3[i][t][0][p][q] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][0][tkj[t]][tkj[k]][tkj[j]];
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        printf("M_tri_18Nn_18Nn: wrong num_derivates\n");
        exit(0);
    }
    
    delete [] tkj;

    // for (int i = 0; i < 1; i++) {
    //     for (int k = 0; k < 18*num_modes*num_scalar; k++) {
    //         for (int j = 0; j < 18*num_modes*num_scalar; j++) {
    //             printf("%0.9f ",M[i][k][j]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
    // }
    // exit(0);

    return M;
}

double ***Mpoisson_tri_18N_18N(double *****MG_tri_6_18_18_18, double ***poisson_tri_18N_10, double ****array_integral4, int num_triangles_in_plane, int num_modes) {
    printf("calculating M2_poisson...\n");
    fflush(stdout); 

    double ***M = new double **[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        M[i] = new double *[18*num_modes];
        for (int k = 0; k < 18*num_modes; k++) {
            M[i][k] = new double [18*num_modes]();
        }
    }

    int *tkj = new int [18*num_modes];
    for (int i = 0; i < 18*num_modes; i++) {
        tkj[i] = i/num_modes;
    }

    #pragma omp parallel for
    for (int i = 0; i < num_triangles_in_plane; i++) {
        for (int t = 0; t < 18*num_modes; t++) {
            for (int k = 0; k < 18*num_modes; k++) {
                for (int j = 0; j < 18*num_modes; j++){
                    // u  v  uu uv vv ζu ζv ζζ ζ
                    // 1  2  3  4  5  6  7  8  9
                    #pragma omp atomic
                    M[i][k][j] += \
                    poisson_tri_18N_10[i][t][0] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][0][tkj[t]][tkj[k]][tkj[j]] + \
                    poisson_tri_18N_10[i][t][1] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][1][tkj[t]][tkj[k]][tkj[j]] + \
                    poisson_tri_18N_10[i][t][2] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][2][tkj[t]][tkj[k]][tkj[j]] + \
                    poisson_tri_18N_10[i][t][3] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][3][tkj[t]][tkj[k]][tkj[j]] + \
                    poisson_tri_18N_10[i][t][4] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][4][tkj[t]][tkj[k]][tkj[j]] + \
                    poisson_tri_18N_10[i][t][5] * array_integral4[t][k][j][0] * MG_tri_6_18_18_18[i][5][tkj[t]][tkj[k]][tkj[j]] + \
                    poisson_tri_18N_10[i][t][9] * array_integral4[t][k][j][1] * MG_tri_6_18_18_18[i][0][tkj[t]][tkj[k]][tkj[j]] + \
                    poisson_tri_18N_10[i][t][6] * array_integral4[t][k][j][1] * MG_tri_6_18_18_18[i][1][tkj[t]][tkj[k]][tkj[j]] + \
                    poisson_tri_18N_10[i][t][7] * array_integral4[t][k][j][1] * MG_tri_6_18_18_18[i][2][tkj[t]][tkj[k]][tkj[j]] + \
                    poisson_tri_18N_10[i][t][8] * array_integral4[t][k][j][2] * MG_tri_6_18_18_18[i][0][tkj[t]][tkj[k]][tkj[j]];
                }
            }
        }
    }
    
    delete [] tkj;

    return M;
}

double **B_tri_18Nn(double ***BG_tri_18_18, double *****object_tri_18N_ope_3_3, int num_triangles_in_plane, int num_modes, int nfp, int num_scalar, double error, int num_derivates) {
    // printf("calculating B2...\n");
    // fflush(stdout); 

    double **B = new double *[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        B[i] = new double [(18*num_modes)*num_scalar]();
    }

    #pragma omp parallel for
    for (int i = 0; i < num_triangles_in_plane; i++) {
        for (int k = 0; k < 18*num_modes; k++) {
            for (int t = 0; t < 18*num_modes; t++) {
                for (int p = 0; p < num_scalar; p++) {
                    #pragma omp atomic
                    B[i][k*num_scalar + p] += object_tri_18N_ope_3_3[i][t][num_derivates][p][0] * integral4(nfp, num_modes, t, k, 0, 0, error) * BG_tri_18_18[i][t/num_modes][k/num_modes];
                }
            }   
        }
    }

    // for (int k = 0; k < 18*num_modes*num_scalar; k++) {
    //     printf("%0.9f\n",B[0][k]);
    // }

    return B;
}

double **Bpoisson_tri_18Nn(double ***BG_tri_18_18, double **object_tri_18N, int num_triangles_in_plane, int num_modes, int nfp, double error) {
    // printf("calculating B2...\n");
    // fflush(stdout); 

    double **B = new double *[num_triangles_in_plane];
    for (int i = 0; i < num_triangles_in_plane; i++) {
        B[i] = new double [18*num_modes]();
    }

    #pragma omp parallel for
    for (int i = 0; i < num_triangles_in_plane; i++) {
        for (int k = 0; k < 18*num_modes; k++) {
            for (int t = 0; t < 18*num_modes; t++) {
                #pragma omp atomic
                B[i][k] += object_tri_18N[i][t] * integral4(nfp, num_modes, t, k, 0, 0, error) * BG_tri_18_18[i][t/num_modes][k/num_modes];
            }   
        }
    }

    return B;
}





double ***solve_surfaces_sparse(double ***B_sec_tri_18, int **points_sequence, int num_triangles, int num_sections) {

    // printf("stiffness matrix assembling...\n");
    // fflush(stdout);

    int num_vertices = points_sequence[0][10];

    VectorXd B(num_vertices*6);
    B.setZero();

    VectorXd X(num_vertices*6);
    X.setZero();

    double ***X_sec_tri_18 = new double **[num_sections];
    for (int i = 0; i < num_sections; i++) {
        X_sec_tri_18[i] = new double *[num_triangles];
        for (int j = 0; j < num_triangles; j++) {
            X_sec_tri_18[i][j] = new double [18]();
        }
    }

    int **index1 = new int *[num_triangles];
    for (int i = 0 ; i < num_triangles; i++) {
        index1[i] = new int [3];
        for (int j = 0; j < 3; j++) {
            index1[i][j] = points_sequence[i][2+3*j] * 6;
        }
    }
    int index2[3] = {0, 6, 12};

    for (int k = 0; k < num_sections; k++) {

        #pragma omp parallel for collapse(3)
        for (int i = 0; i < num_triangles; i++) {
            for (int j = 0; j < 3; j++) {
                for (int s = 0; s < 6; s++) {
                    #pragma omp atomic
                    B(index1[i][j] + s) += B_sec_tri_18[k][i][index2[j] + s];
                }
            }
        }

        X = solver_LU.solve(B);

        #pragma omp parallel for collapse(3)
        for (int i = 0; i < num_triangles; i++) {
            for (int j = 0; j < 3; j++) {
                for (int s = 0; s < 6; s++) {
                    X_sec_tri_18[k][i][index2[j] + s] = X[index1[i][j] + s];
                }
            }
        }

        B.setZero();

    }

    for (int i = 0 ; i < num_triangles; i++) {
        delete [] index1[i];
    }
    delete [] index1;
    
    return X_sec_tri_18;
}

double ****solve_surfaces_sparse_num(double ****obe_nu_sec_tri_18_num, double ***temp_sec_tri_18, double ***B_sec_tri_18, int **points_sequence, int num_sections, int num_triangles, int num_d) {

    double ****obe_sec_tri_18_num = new double ***[num_sections];
    for (int i = 0; i < num_sections; i++) {
        obe_sec_tri_18_num[i] = new double **[num_triangles];
        for (int j = 0; j < num_triangles; j++) {
            obe_sec_tri_18_num[i][j] = new double *[18];
            for (int k = 0; k < 18; k++) {
                obe_sec_tri_18_num[i][j][k] = new double [num_d]();
            }
        }
    }
    
    for (int i = 0; i < num_d; i++) {
        #pragma omp parallel for
        for (int l = 0; l < num_sections; l++) {
            for (int s = 0; s < num_triangles; s++) {
                for (int t = 0; t < 18; t++) {
                    B_sec_tri_18[l][s][t] = obe_nu_sec_tri_18_num[l][s][t][i];
                }
            }
        }

        temp_sec_tri_18 = solve_surfaces_sparse(B_sec_tri_18, points_sequence, num_triangles, num_sections);

        #pragma omp parallel for
        for (int l = 0; l < num_sections; l++) {
            for (int s = 0; s < num_triangles; s++) {
                for (int t = 0; t < 18; t++) {
                    obe_sec_tri_18_num[l][s][t][i] = temp_sec_tri_18[l][s][t];
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

    return obe_sec_tri_18_num;
}

double ***solve_magnetic_field_sparse_direct(double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangels_in_plane, int num_scalar, double error, time_t time_s) {
    printf("stiffness matrix assembling...\n");
    fflush(stdout);
    int num_edges = points_sequence[0][9];
    int num_vertices = points_sequence[0][10];
    // printf("num_edges = %d\n", num_edges);

    int num_modes_scalar = num_modes*num_scalar;
    int num_total_DoFs = num_vertices*6*num_modes_scalar; 

    int num_free_DoFs = (num_vertices*6-num_edges)*num_modes_scalar;
    int num_positive_DoFs = (num_vertices-num_edges)*6*num_modes_scalar;
    // int num_negative_DoFs = num_edges*5*num_modes_scalar;
    int num_fixed_DoFs = num_edges*num_modes_scalar;
    
    
    double **X = new double *[num_scalar];
    for (int i = 0; i < num_scalar; i++) {
        X[i] = new double [num_vertices*6*num_modes]();
    }

    SparseMatrix<double> A11(num_free_DoFs, num_free_DoFs); 
    SparseMatrix<double> A12(num_free_DoFs, num_fixed_DoFs); 
    A11.reserve(num_vertices*6*num_modes_scalar*6*num_modes_scalar - num_edges*num_modes_scalar*num_modes_scalar); // (num_vertices-num_edges)*6*num_modes
    A12.reserve((num_vertices-num_edges)*6*num_modes_scalar*6*num_modes_scalar);
    A11.setZero();
    A12.setZero();

    VectorXd X1(num_free_DoFs);
    VectorXd X2(num_fixed_DoFs);
    VectorXd M(num_free_DoFs);
    X1.setZero();
    X2.setZero();
    M.setZero();

    VectorXd B1(num_free_DoFs);
    B1.setZero();

    for (int i = 0; i < num_fixed_DoFs; i++) {
        X2(i) = boundary_DoFs[i%num_scalar][i/num_scalar];
    }

    // Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver1;
    // Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver1;
    // Eigen::LeastSquaresConjugateGradient<SparseMatrix<double>>> solver1;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver2;
    // Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver2;
    
    int flag_j, flag_k, index_J, index_K, index_j1, index_k1, index_k2;

    typedef Eigen::Triplet<double> T;
    vector<T> coeff1, coeff2;

    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            flag_j = points_sequence[i][0+3*j];
            index_J = 6*num_modes_scalar*j;
            for (int k = 0; k < 3; k++) {
                flag_k = points_sequence[i][0+3*k];
                index_K = 6*num_modes_scalar*k;

                if (flag_j == 0) {
                    index_j1 = points_sequence[i][1+3*j]*6*num_modes_scalar;
                    if (flag_k == 0) {
                        index_k1 = points_sequence[i][1+3*k]*6*num_modes_scalar;
                        for (int s = 0; s < 6*num_modes_scalar; s++) {
                            for (int t = 0; t < 6*num_modes_scalar; t++) {
                                coeff1.push_back( T(index_j1 + s, index_k1 + t, M_tri_18Nn_18Nn[i][index_J + s][index_K + t]) );
                            }
                        }
                    }
                    else if (flag_k == 1) {
                        index_k1 = num_positive_DoFs + points_sequence[i][1+3*k]*5*num_modes_scalar;
                        index_k2 = points_sequence[i][1+3*k]*num_modes_scalar;
                        for (int s = 0; s < 6*num_modes_scalar; s++) {
                            for (int t = 0; t < 5*num_modes_scalar; t++) {
                                coeff1.push_back( T(index_j1 + s, index_k1 + t, M_tri_18Nn_18Nn[i][index_J + s][index_K + num_modes_scalar + t]) );
                            }
                        }
                        for (int s = 0; s < 6*num_modes_scalar; s++) {
                            for (int t = 0; t < num_modes_scalar; t++) {
                                coeff2.push_back( T(index_j1 + s, index_k2 + t, M_tri_18Nn_18Nn[i][index_J + s][index_K + t]) );
                            }
                        }
                    }
                }
                else if (flag_j == 1) {
                    index_j1 = num_positive_DoFs + points_sequence[i][1+3*j]*5*num_modes_scalar;
                    if (flag_k == 0) {
                        index_k1 = points_sequence[i][1+3*k]*6*num_modes_scalar;
                        for (int s = 0; s < 5*num_modes_scalar; s++) {
                            for (int t = 0; t < 6*num_modes_scalar; t++) {
                                coeff1.push_back( T(index_j1 + s, index_k1 + t, M_tri_18Nn_18Nn[i][index_J + num_modes_scalar + s][index_K + t]) );
                            }
                        }
                    }
                    else if (flag_k == 1) {
                        index_k1 = num_positive_DoFs + points_sequence[i][1+3*k]*5*num_modes_scalar;
                        index_k2 = points_sequence[i][1+3*k]*num_modes_scalar;
                        for (int s = 0; s < 5*num_modes_scalar; s++) {
                            for (int t = 0; t < 5*num_modes_scalar; t++) {
                                coeff1.push_back( T(index_j1 + s, index_k1 + t, M_tri_18Nn_18Nn[i][index_J + num_modes_scalar + s][index_K + num_modes_scalar+ t]) );
                            }
                        }
                        for (int s = 0; s < 5*num_modes_scalar; s++) {
                            for (int t = 0; t < num_modes_scalar; t++) {
                                coeff2.push_back( T(index_j1 + s, index_k2 + t, M_tri_18Nn_18Nn[i][index_J + num_modes_scalar + s][index_K + t]) );
                            }
                        }
                    }
                } 

            }
        }
    }

    A11.setFromTriplets(coeff1.begin(), coeff1.end());
    A12.setFromTriplets(coeff2.begin(), coeff2.end());

    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            flag_j = points_sequence[i][0+3*j];
            index_J = 6*num_modes_scalar*j;
            if (flag_j == 0) {
                index_j1 = points_sequence[i][1+3*j]*6*num_modes_scalar;
                for (int s = 0; s < 6*num_modes_scalar; s++) {
                    B1(index_j1 + s) += B_tri_18Nn[i][index_J + s];
                }
            }
            else if (flag_j == 1) {
                index_j1 = num_positive_DoFs + points_sequence[i][1+3*j]*5*num_modes_scalar;
                for (int s = 0; s < 5*num_modes_scalar; s++) {
                    B1(index_j1 + s) += B_tri_18Nn[i][index_J + num_modes_scalar + s];
                }
            }
        }
    }

    M = B1 - A12 * X2;

    A12.resize(1,1);
    A12.data().squeeze();
    B1.resize(1);
    
    timing(time_s);

    printf("solving linear equation...\n");
    fflush(stdout);    

    // printf(" preconditioning for solver...\n");
    // fflush(stdout);
    // A11.makeCompressed();
    // solver1.compute(A11);
    // if (solver1.info() != Eigen::Success) {
    //     printf("solve_beltrami_sparse: BiCGSTAB failed...\n");
    //     exit(0);
    // }
    // timing(time_s);

    // solver1.setMaxIterations(num_edges*10*num_scalar);
    // solver1.setTolerance(error);
    // X1 = solver1.solve(M);
    // // X1 = solver1.solveWithGuess(M, X1_guess);
    // std::cout << " #iterations:     " << solver1.iterations() << std::endl;
    // std::cout << " #estimated error: " << solver1.error() << std::endl;    

    printf(" preconditioning...\n");
    fflush(stdout); 
    A11.makeCompressed();   
    solver2.compute(A11);
    if (solver2.info() != Eigen::Success) {
        printf("solve_beltrami_sparse: spraseLU failed...\n");
        exit(0);
    }
    timing(time_s);

    X1 = solver2.solve(M);
    std::cout << " #estimated error: " << (A11 * X1 - M).norm() / M.norm() << std::endl;

    printf(" renumbering...\n");
    fflush(stdout);

    for (int i = 0; i < num_free_DoFs; i++) {
        if (i < num_positive_DoFs) {
            X[i%num_scalar][i/num_scalar] = X1(i);
        }
        else {
            X[i%num_scalar][(i + ((i-num_positive_DoFs) / (5*num_modes_scalar) + 1) * num_modes_scalar) / num_scalar] = X1(i);
        }
    }
    for (int i = 0; i < num_fixed_DoFs; i++) {
        X[i%num_scalar][(i + num_positive_DoFs + (i/num_modes_scalar) * (5*num_modes_scalar)) / num_scalar] = X2(i);
    }

    // for (int j = 0; j < num_scalar; j++) {
    //     for (int i = 0; i < num_total_DoFs/num_scalar; i++) {
    //         printf("%0.9f  ", X[j][i]);
    //         if ((i+1)%(num_modes) == 0) {
    //             printf("        ");
    //         }
    //         if ((i+1)%(6*num_modes) == 0) {
    //             printf("\n");
    //         }
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // // exit(0);

    int kk = 0;
    double ***F_sca_tri_18N = new double **[num_scalar];
    for (int i = 0; i < num_scalar; i++) {
        F_sca_tri_18N[i] = new double *[num_triangels_in_plane];
        for (int j = 0; j < num_triangels_in_plane; j++) {
            F_sca_tri_18N[i][j] = new double [18*num_modes] ();
            for (int k = 0; k < 18*num_modes; k++) {
                kk = k/(6*num_modes);
                if (points_sequence[j][0+3*kk] == 0) {
                    F_sca_tri_18N[i][j][k] = X[i][points_sequence[j][1+3*kk]*6*num_modes + k - kk*6*num_modes];
                }
                else if (points_sequence[j][0+3*kk] == 1) {
                    F_sca_tri_18N[i][j][k] = X[i][num_positive_DoFs/num_scalar + points_sequence[j][1+3*kk]*6*num_modes + k - kk*6*num_modes];
                }
                else {
                    printf("solve_1d: error in F_tri_18N...\n");
                    exit(0);
                }
            }
        }
    } 

    delete [] X;
    
    return F_sca_tri_18N;   
}

double ***solve_magnetic_field_sparse_guess(double ***guess_sca_tri_18N, double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangels_in_plane, int num_scalar, double error, time_t time_s) {
    printf("stiffness matrix assembling...\n");
    fflush(stdout);
    int num_edges = points_sequence[0][9];
    int num_vertices = points_sequence[0][10];
    // printf("num_edges = %d\n", num_edges);

    int num_modes_scalar = num_modes*num_scalar;
    int num_total_DoFs = num_vertices*6*num_modes_scalar; 

    int num_free_DoFs = (num_vertices*6-num_edges)*num_modes_scalar;
    int num_positive_DoFs = (num_vertices-num_edges)*6*num_modes_scalar;
    // int num_negative_DoFs = num_edges*5*num_modes_scalar;
    int num_fixed_DoFs = num_edges*num_modes_scalar;
    
    
    double **X = new double *[num_scalar];
    for (int i = 0; i < num_scalar; i++) {
        X[i] = new double [num_vertices*6*num_modes]();
    }

    SparseMatrix<double> A11(num_free_DoFs, num_free_DoFs); 
    SparseMatrix<double> A12(num_free_DoFs, num_fixed_DoFs); 
    A11.reserve(num_vertices*6*num_modes_scalar*6*num_modes_scalar - num_edges*num_modes_scalar*num_modes_scalar); // (num_vertices-num_edges)*6*num_modes
    A12.reserve((num_vertices-num_edges)*6*num_modes_scalar*6*num_modes_scalar);
    A11.setZero();
    A12.setZero();

    VectorXd X1(num_free_DoFs);
    VectorXd X1_guess(num_free_DoFs);
    VectorXd X2(num_fixed_DoFs);
    VectorXd M(num_free_DoFs);
    X1.setZero();
    X1_guess.setZero();
    X2.setZero();
    M.setZero();

    VectorXd B1(num_free_DoFs);
    B1.setZero();

    for (int i = 0; i < num_fixed_DoFs; i++) {
        X2(i) = boundary_DoFs[i%num_scalar][i/num_scalar];
        // X2(i) = 0.0;
    }

    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            if (points_sequence[i][0+3*j] == 0) {
                for (int s = 0; s < 6*num_modes_scalar; s++) {
                    X1_guess(points_sequence[i][1+3*j]*6*num_modes_scalar + s) = guess_sca_tri_18N[s%num_scalar][i][j * (6*num_modes) + s/num_scalar];
                }
            }
            else if (points_sequence[i][0+3*j] == 1) {
                for (int s = 0; s < 5*num_modes_scalar; s++) {
                    X1_guess(num_positive_DoFs + points_sequence[i][1+3*j]*5*num_modes_scalar + s) = guess_sca_tri_18N[s%num_scalar][i][j * (6*num_modes) + num_modes + s/num_scalar];
                }
            }
            else {
                printf("solve_magnetic_field_sparse_direct_guess: incorrect sequence\n");
                fflush(stdout);
                exit(0);
            }
        }
    }

    // Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver1;
    // Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver1;
    Eigen::LeastSquaresConjugateGradient<SparseMatrix<double>> solver1;
    // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver2;
    // Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver2;
    
    int flag_j, flag_k, index_J, index_K, index_j1, index_k1, index_k2;

    typedef Eigen::Triplet<double> T;
    vector<T> coeff1, coeff2;

    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            flag_j = points_sequence[i][0+3*j];
            index_J = 6*num_modes_scalar*j;
            for (int k = 0; k < 3; k++) {
                flag_k = points_sequence[i][0+3*k];
                index_K = 6*num_modes_scalar*k;

                if (flag_j == 0) {
                    index_j1 = points_sequence[i][1+3*j]*6*num_modes_scalar;
                    if (flag_k == 0) {
                        index_k1 = points_sequence[i][1+3*k]*6*num_modes_scalar;
                        for (int s = 0; s < 6*num_modes_scalar; s++) {
                            for (int t = 0; t < 6*num_modes_scalar; t++) {
                                coeff1.push_back( T(index_j1 + s, index_k1 + t, M_tri_18Nn_18Nn[i][index_J + s][index_K + t]) );
                            }
                        }
                    }
                    else if (flag_k == 1) {
                        index_k1 = num_positive_DoFs + points_sequence[i][1+3*k]*5*num_modes_scalar;
                        index_k2 = points_sequence[i][1+3*k]*num_modes_scalar;
                        for (int s = 0; s < 6*num_modes_scalar; s++) {
                            for (int t = 0; t < 5*num_modes_scalar; t++) {
                                coeff1.push_back( T(index_j1 + s, index_k1 + t, M_tri_18Nn_18Nn[i][index_J + s][index_K + num_modes_scalar + t]) );
                            }
                        }
                        for (int s = 0; s < 6*num_modes_scalar; s++) {
                            for (int t = 0; t < num_modes_scalar; t++) {
                                coeff2.push_back( T(index_j1 + s, index_k2 + t, M_tri_18Nn_18Nn[i][index_J + s][index_K + t]) );
                            }
                        }
                    }
                }
                else if (flag_j == 1) {
                    index_j1 = num_positive_DoFs + points_sequence[i][1+3*j]*5*num_modes_scalar;
                    if (flag_k == 0) {
                        index_k1 = points_sequence[i][1+3*k]*6*num_modes_scalar;
                        for (int s = 0; s < 5*num_modes_scalar; s++) {
                            for (int t = 0; t < 6*num_modes_scalar; t++) {
                                coeff1.push_back( T(index_j1 + s, index_k1 + t, M_tri_18Nn_18Nn[i][index_J + num_modes_scalar + s][index_K + t]) );
                            }
                        }
                    }
                    else if (flag_k == 1) {
                        index_k1 = num_positive_DoFs + points_sequence[i][1+3*k]*5*num_modes_scalar;
                        index_k2 = points_sequence[i][1+3*k]*num_modes_scalar;
                        for (int s = 0; s < 5*num_modes_scalar; s++) {
                            for (int t = 0; t < 5*num_modes_scalar; t++) {
                                coeff1.push_back( T(index_j1 + s, index_k1 + t, M_tri_18Nn_18Nn[i][index_J + num_modes_scalar + s][index_K + num_modes_scalar+ t]) );
                            }
                        }
                        for (int s = 0; s < 5*num_modes_scalar; s++) {
                            for (int t = 0; t < num_modes_scalar; t++) {
                                coeff2.push_back( T(index_j1 + s, index_k2 + t, M_tri_18Nn_18Nn[i][index_J + num_modes_scalar + s][index_K + t]) );
                            }
                        }
                    }
                } 

            }
        }
    }

    A11.setFromTriplets(coeff1.begin(), coeff1.end());
    A12.setFromTriplets(coeff2.begin(), coeff2.end());

    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            flag_j = points_sequence[i][0+3*j];
            index_J = 6*num_modes_scalar*j;
            if (flag_j == 0) {
                index_j1 = points_sequence[i][1+3*j]*6*num_modes_scalar;
                for (int s = 0; s < 6*num_modes_scalar; s++) {
                    B1(index_j1 + s) += B_tri_18Nn[i][index_J + s];
                }
            }
            else if (flag_j == 1) {
                index_j1 = num_positive_DoFs + points_sequence[i][1+3*j]*5*num_modes_scalar;
                for (int s = 0; s < 5*num_modes_scalar; s++) {
                    B1(index_j1 + s) += B_tri_18Nn[i][index_J + num_modes_scalar + s];
                }
            }
        }
    }

    M = B1 - A12 * X2;

    A12.resize(1,1);
    A12.data().squeeze();
    B1.resize(1);
    
    timing(time_s);

    printf("solving linear equation...\n");
    fflush(stdout);    

    printf(" preconditioning for solver...\n");
    fflush(stdout);
    A11.makeCompressed();
    solver1.compute(A11);
    if (solver1.info() != Eigen::Success) {
        printf("solve_beltrami_sparse: BiCGSTAB failed...\n");
        exit(0);
    }
    timing(time_s);

    // solver1.setMaxIterations(num_edges*10*num_scalar);
    solver1.setMaxIterations(1999);
    solver1.setTolerance(error);
    // X1 = solver1.solve(M);
    X1 = solver1.solveWithGuess(M, X1_guess);
    std::cout << " #iterations:     " << solver1.iterations() << std::endl;
    std::cout << " #estimated error: " << solver1.error() << std::endl;    

    // printf(" preconditioning...\n");
    // fflush(stdout); 
    // A11.makeCompressed();   
    // solver2.compute(A11);
    // if (solver2.info() != Eigen::Success) {
    //     printf("solve_beltrami_sparse: spraseLU failed...\n");
    //     exit(0);
    // }
    // timing(time_s);

    // X1 = solver2.solve(M);
    // std::cout << " #estimated error: " << (A11 * X1 - M).norm() / M.norm() << std::endl;

    // printf(" renumbering...\n");
    // fflush(stdout);

    for (int i = 0; i < num_free_DoFs; i++) {
        if (i < num_positive_DoFs) {
            X[i%num_scalar][i/num_scalar] = X1(i);
        }
        else {
            X[i%num_scalar][(i + ((i-num_positive_DoFs) / (5*num_modes_scalar) + 1) * num_modes_scalar) / num_scalar] = X1(i);
        }
    }
    for (int i = 0; i < num_fixed_DoFs; i++) {
        X[i%num_scalar][(i + num_positive_DoFs + (i/num_modes_scalar) * (5*num_modes_scalar)) / num_scalar] = X2(i);
    }

    // for (int j = 0; j < num_scalar; j++) {
    //     for (int i = 0; i < num_total_DoFs/num_scalar; i++) {
    //         printf("%0.9f  ", X[j][i]);
    //         if ((i+1)%(num_modes) == 0) {
    //             printf("        ");
    //         }
    //         if ((i+1)%(6*num_modes) == 0) {
    //             printf("\n");
    //         }
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // // exit(0);

    int kk = 0;
    double ***F_sca_tri_18N = new double **[num_scalar];
    for (int i = 0; i < num_scalar; i++) {
        F_sca_tri_18N[i] = new double *[num_triangels_in_plane];
        for (int j = 0; j < num_triangels_in_plane; j++) {
            F_sca_tri_18N[i][j] = new double [18*num_modes] ();
            for (int k = 0; k < 18*num_modes; k++) {
                kk = k/(6*num_modes);
                if (points_sequence[j][0+3*kk] == 0) {
                    F_sca_tri_18N[i][j][k] = X[i][points_sequence[j][1+3*kk]*6*num_modes + k - kk*6*num_modes];
                }
                else if (points_sequence[j][0+3*kk] == 1) {
                    F_sca_tri_18N[i][j][k] = X[i][num_positive_DoFs/num_scalar + points_sequence[j][1+3*kk]*6*num_modes + k - kk*6*num_modes];
                }
                else {
                    printf("solve_1d: error in F_tri_18N...\n");
                    exit(0);
                }
            }
        }
    } 

    delete [] X;
    
    return F_sca_tri_18N;   
}

double ***solve_magnetic_field_sparse_lagrange_multiplier(double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangels_in_plane, int num_scalar, int flag, double error, time_t time_s) {

    printf("stiffness matrix assembling...\n");
    fflush(stdout);

    int num_edges = points_sequence[0][9];
    int num_vertices = points_sequence[0][10];
    // printf("num_edges = %d\n", num_edges);

    int num_modes_scalar = num_modes*num_scalar;
    int num_total_DoFs = num_vertices*6*num_modes_scalar;
    int num_boundary_DoFs = num_edges*num_modes_scalar;
    if (flag == 1) {
        if (num_scalar != 1) {
            printf("num_scalar should equal to 1 for divergence cleaning...\n");
            fflush(stdout);
            exit(0);
        }
        num_boundary_DoFs = num_edges*(3*num_modes - 1);
    }

    if (flag == 2) {
        num_boundary_DoFs = 0;
    }

    int num_extended_DoFs = num_total_DoFs + num_boundary_DoFs;
    if (flag == 1) {
        num_extended_DoFs += 1;
    }
    
    SparseMatrix<double> A(num_extended_DoFs, num_extended_DoFs);
    A.setZero();

    VectorXd B(num_extended_DoFs);
    B.setZero();

    VectorXd X(num_extended_DoFs);
    VectorXd X_upper(num_total_DoFs);
    X.setZero();
    X_upper.setZero();

    double **X_re = new double *[num_scalar];
    for (int i = 0; i < num_scalar; i++) {
        X_re[i] = new double [num_vertices*6*num_modes]();
    }

    typedef Eigen::Triplet<double> T;
    vector<T> coeff_A, coeff_A_asy;

    // Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver1;
    Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver1;
    // Eigen::LeastSquaresConjugateGradient<SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver1;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver2;
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver3;

    int index_J, index_K, index_j1, index_j2, index_k1;

    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            index_J = 6*num_modes_scalar*j;
            index_j1 = points_sequence[i][2+3*j]*6*num_modes_scalar;
            for (int k = 0; k < 3; k++) {
                index_K = 6*num_modes_scalar*k;
                index_k1 = points_sequence[i][2+3*k]*6*num_modes_scalar;
                for (int s = 0; s < 6*num_modes_scalar; s++) {
                    for (int t = 0; t < 6*num_modes_scalar; t++) {
                        coeff_A.push_back( T(index_j1 + s, index_k1 + t, M_tri_18Nn_18Nn[i][index_J + s][index_K + t]) );
                    }
                }
            }
        }
    }

    vector<int> vec2;
    int flag2 = 0;
    int index1, index2, index3;
    int point_on_boundary;
    for (int i = 0; i < num_triangels_in_plane; i++) {
        // dirichlet boundary condition 
        // or neumann boundary condition for divergence cleaning
        // or skip for projection2
        if (flag == 2) {
            continue;
        }

        for (int j = 0; j < 3; j++) {    
            if (points_sequence[i][0+3*j] == 1) {

                point_on_boundary = points_sequence[i][2+3*j];

                for (int k = 0; k < vec2.size(); k++) {
                    if (vec2[k] == points_sequence[i][1+3*j]) {
                        flag2 = 1;
                        break;
                    }
                }
                if (flag2 == 1) {
                    flag2 = 0;
                    continue;
                }
                else {
                    vec2.push_back(points_sequence[i][1+3*j]);
                }

                if (flag == 0) {
                    index2 = points_sequence[i][1+3*j]*num_modes_scalar;
                    index3 = points_sequence[i][2+3*j]*6*num_modes_scalar;
                    for (int s = 0; s < num_modes; s++) {
                        for (int t = 0; t < num_scalar; t++) {
                            coeff_A.push_back( T(num_total_DoFs + index2 + num_scalar*s + t, index3 + num_scalar*s + t, 1.0) );
                            coeff_A.push_back( T(index3 + num_scalar*s + t, num_total_DoFs + index2 + num_scalar*s + t, 1.0) );
                        }
                    }
                }
                else if (flag == 1) {
                    index2 = points_sequence[i][1+3*j]*(3*num_modes - 1);
                    index3 = points_sequence[i][2+3*j]*(6*num_modes) + 1;
                    for (int s = 0; s < 3*num_modes - 1; s++) {    
                        coeff_A.push_back( T(num_total_DoFs + index2 + s, index3 + s, 1.0) );
                        coeff_A.push_back( T(index3 + s, num_total_DoFs + index2 + s, 1.0) );
                    }
                }
                else {
                    printf("solve_magnetic_field_sparse_lagrange_multiplier: error in boundary condition...\n");
                    exit(0);
                }
            }
        }
    }

    if (flag == 1) {
        coeff_A.push_back( T(num_total_DoFs + num_boundary_DoFs, point_on_boundary*(6*num_modes), 1.0) );
        coeff_A.push_back( T(point_on_boundary*(6*num_modes), num_total_DoFs + num_boundary_DoFs, 1.0) );
        B(num_total_DoFs + num_boundary_DoFs) = 1.0;
    }


    A.setFromTriplets(coeff_A.begin(), coeff_A.end());
    A.makeCompressed();
    vector<T>().swap(coeff_A);

    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            index1 = 6*num_modes_scalar*j;
            index2 = points_sequence[i][2+3*j]*6*num_modes_scalar;
            for (int s = 0; s < 6*num_modes_scalar; s++) {
                B(index2 + s) += B_tri_18Nn[i][index1 + s];
            }

            if (points_sequence[i][0+3*j] == 1 && flag == 0) {
                index2 = points_sequence[i][1+3*j]*num_modes_scalar;
                index3 = points_sequence[i][1+3*j]*num_modes;
                for (int s = 0; s < num_modes; s++) {
                    for (int t = 0; t < num_scalar; t++) { // num_edges*num_modes
                        B(num_total_DoFs + index2 + num_scalar*s + t) = boundary_DoFs[t][index3 + s];
                    }
                }
            }

        }
    }

    timing(time_s);

    printf("solving linear equation...\n");
    fflush(stdout);    

    printf(" preconditioning...\n");
    solver2.compute(A);
    if (solver2.info() != Eigen::Success) {
        printf("solve_magnetic_field_sparse_lagrange_multiplier: decomposition failed...\n");
        exit(0);
    }
    fflush(stdout); 

    X = solver2.solve(B);
    printf("estimated error = %0.9f\n", abs((A * X).norm() - B.norm()) / (B.norm()));
    timing(time_s);

    // printf(" preconditioning...\n");
    // solver1.compute(A);
    // if (solver1.info() != Eigen::Success) {
    //     printf("solve_magnetic_field_sparse_lagrange_multiplier: decomposition failed...\n");
    //     exit(0);
    // }
    // fflush(stdout); 
    // solver1.setMaxIterations(1999);
    // // X = solver1.solveWithGuess(B, X.eval());
    // X = solver1.solve(B);
    // printf("estimated error = %0.9f, iterations = %d\n", solver1.error(), solver1.iterations());
    // timing(time_s);

    for (int i = 0; i < num_total_DoFs; i++) {
        X_upper(i) = X(i);
    }


    printf(" renumbering...\n");
    fflush(stdout);
    for (int i = 0; i < num_total_DoFs; i++) {
        X_re[i%num_scalar][i/num_scalar] = X_upper(i);
    }

    // for (int j = 0; j < num_scalar; j++) {
    //     for (int i = 0; i < num_total_DoFs/num_scalar; i++) {
    //         printf("%0.9f  ", X_re[j][i]);
    //         if ((i+1)%(num_modes) == 0) {
    //             printf("        ");
    //         }
    //         if ((i+1)%(6*num_modes) == 0) {
    //             printf("\n");
    //         }
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // fflush(stdout);
    // // exit(0);

    int kk = 0;
    double ***F_sca_tri_18N = new double **[num_scalar];
    for (int i = 0; i < num_scalar; i++) {
        F_sca_tri_18N[i] = new double *[num_triangels_in_plane];
        for (int j = 0; j < num_triangels_in_plane; j++) {
            F_sca_tri_18N[i][j] = new double [18*num_modes] ();
            for (int k = 0; k < 18*num_modes; k++) {
                kk = k/(6*num_modes);
                F_sca_tri_18N[i][j][k] = X_re[i][points_sequence[j][2+3*kk]*6*num_modes + k - kk*6*num_modes];
            }
        }
    } 

    delete [] X_re;    

    return F_sca_tri_18N;

}

double ***solve_magnetic_field_sparse_lagrange_multiplier_div(double**** divergence_tri_3_var_modes, double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangels_in_plane, int nfp_ft, int num_scalar, int flag, double error, time_t time_s) {

    printf("stiffness matrix assembling...\n");
    fflush(stdout);

    int num_variables = 9;
    int num_edges = points_sequence[0][9];
    int num_vertices = points_sequence[0][10];
    // int num_modes_divergence = 2*num_modes-1;
    int num_modes_divergence = 0;
    int num_modes_nromal = (2*num_modes-1);
    num_modes_nromal = 1;
    // printf("num_edges = %d\n", num_edges);

    int num_modes_scalar = num_modes*num_scalar;
    int num_total_DoFs = num_vertices * (6*num_modes_scalar);
    // int num_divergence_DoFs = 0;
    int num_divergence_DoFs = num_vertices * num_modes_divergence;
    int num_boundary_DoFs = num_edges * num_modes_scalar;
    
    int num_extended_DoFs = num_total_DoFs + num_boundary_DoFs;
    
    SparseMatrix<double> A(num_extended_DoFs + num_divergence_DoFs, num_extended_DoFs + num_divergence_DoFs);
    A.setZero();

    VectorXd B(num_extended_DoFs + num_divergence_DoFs);
    B.setZero();

    VectorXd X(num_extended_DoFs + num_divergence_DoFs);
    VectorXd X_upper(num_total_DoFs);
    X.setZero();
    X_upper.setZero();

    double **X_re = new double *[num_scalar];
    for (int i = 0; i < num_scalar; i++) {
        X_re[i] = new double [num_vertices*6*num_modes]();
    }

    typedef Eigen::Triplet<double> T;
    vector<T> coeff_A;

    // Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver2;
    // Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver1;
    // Eigen::LeastSquaresConjugateGradient<SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver1;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver2;
    // Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver2;

    int index_J, index_K, index_j1, index_j2, index_k1;

    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            index_J = 6*num_modes_scalar*j;
            index_j1 = points_sequence[i][2+3*j]*6*num_modes_scalar;
            for (int k = 0; k < 3; k++) {
                index_K = 6*num_modes_scalar*k;
                index_k1 = points_sequence[i][2+3*k]*6*num_modes_scalar;
                for (int s = 0; s < 6*num_modes_scalar; s++) {
                    for (int t = 0; t < 6*num_modes_scalar; t++) {
                        coeff_A.push_back( T(index_j1 + s, index_k1 + t, M_tri_18Nn_18Nn[i][index_J + s][index_K + t]) );
                    }
                }
            }
        }
    }

    int *p1;
    double *p2;
    vector<int> vec2;
    int flag2 = 0;
    int index1, index2, index3, index4;
    // dirichlet boundary condition 
    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {    
            if (points_sequence[i][0+3*j] == 1) {

                for (int k = 0; k < vec2.size(); k++) {
                    if (vec2[k] == points_sequence[i][1+3*j]) {
                        flag2 = 1;
                        break;
                    }
                }
                if (flag2 == 1) {
                    flag2 = 0;
                    continue;
                }
                else {
                    vec2.push_back(points_sequence[i][1+3*j]);
                }

                index2 = points_sequence[i][1+3*j]*num_modes_scalar;
                index3 = points_sequence[i][2+3*j]*6*num_modes_scalar;
                for (int s = 0; s < num_modes; s++) {
                    for (int t = 0; t < num_scalar; t++) {
                        coeff_A.push_back( T(num_total_DoFs + index2 + num_scalar*s + t, index3 + num_scalar*s + t, 1.0) );
                        coeff_A.push_back( T(index3 + num_scalar*s + t, num_total_DoFs + index2 + num_scalar*s + t, 1.0) );
                    }
                }
                
            }

        }
    }

    
    std::function<int*(int, int)> position_s = [&] (int s, int t) {
        int *pos = new int [2]();

        int mode_s;
        int mode_t;
        if (s % 2 == 1) {
            mode_s = (s + 1) / 2;
        }
        else {
            mode_s = s / 2;
        }
        if (t % 2 == 1) {
            mode_t = (t + 1) / 2;
        }
        else {
            mode_t = t / 2;
        }

        if (s == 0) { // 0_X
            pos[0] = t;
            pos[1] = t;
        }
        else if (t == 0) { // X_0
            pos[0] = s;
            pos[1] = s;
        }
        else if (s % 2 == 1) {
            if (t % 2 == 1) { // cos_cos
                pos[0] = (mode_s + mode_t) * 2 - 1;
                if (mode_s == mode_t) {
                    pos[1] = 0;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2 - 1;
                }
            }
            else if (t % 2 == 0) { // cos_sin
                pos[0] = (mode_s + mode_t) * 2;
                if (mode_s == mode_t) {
                    pos[1] = -999;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2;
                }
            }
        }
        else if (s % 2 == 0) {
            if (t % 2 == 1) { // sin_cos
                pos[0] = (mode_s + mode_t) * 2;
                if (mode_s == mode_t) {
                    pos[1] = -999;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2;
                }
            }
            else if (t % 2 == 0) { // sin_sin
                pos[0] = (mode_s + mode_t) * 2 - 1;
                if (mode_s == mode_t) {
                    pos[1] = 0;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2 - 1;
                }
            }
        }

        return pos;
    };

    std::function<double*(int, int, int)> coefficients_v = [&] (int s, int t, int nfp_ft) {
        double *coeff = new double [2]();

        double sign = -999999;
        if ((s - t) != 0) {
            sign = (s - t) / abs(s - t);
        }

        if (s == 0 || t == 0) { // 0_X X_0
            coeff[0] = 0.5;
            coeff[1] = 0.5;
        }
        else if (s % 2 == 1) {
            if (t % 2 == 1) { // cos_cos
                coeff[0] =  0.5;
                coeff[1] =  0.5;
            }
            else if (t % 2 == 0) { // cos_sin
                coeff[0] =  0.5;
                coeff[1] = -0.5 * sign;
            }
        }
        else if (s % 2 == 0) {
            if (t % 2 == 1) { // sin_cos
                coeff[0] =  0.5;
                coeff[1] =  0.5 * sign;
            }
            else if (t % 2 == 0) { // sin_sin
                coeff[0] = -0.5;
                coeff[1] =  0.5;
            }
        }

        return coeff;
    };

    std::function<int*(int, int)> position_s_ex = [&] (int s, int t) {
        int *pos = new int [2]();

        int mode_s;
        int mode_t;
        if (s % 2 == 1) {
            mode_s = (s + 1) / 2;
        }
        else {
            mode_s = s / 2;
        }
        if (t % 2 == 1) {
            mode_t = (t + 1) / 2;
        }
        else {
            mode_t = t / 2;
        }
        
        if (s == 0) { // 0_X
            pos[0] = -999;
            pos[1] = -999;
        }
        else if (t == 0) { // X_0
            if (s % 2 == 1) {
                pos[0] = s + 1;
                pos[1] = s + 1;
            }
            else if (s % 2 == 0) {
                pos[0] = s - 1;
                pos[1] = s - 1;
            }
        }
        else if (s % 2 == 1) {
            if (t % 2 == 1) { // sin_cos
                pos[0] = (mode_s + mode_t) * 2;
                if (mode_s == mode_t) {
                    pos[1] = -999;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2;
                }
            }
            else if (t % 2 == 0) { // sin_sin
                pos[0] = (mode_s + mode_t) * 2 - 1;
                if (mode_s == mode_t) {
                    pos[1] = 0;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2 - 1;
                }
            }
        }
        else if (s % 2 == 0) {
            if (t % 2 == 1) { // cos_cos
                pos[0] = (mode_s + mode_t) * 2 - 1;
                if (mode_s == mode_t) {
                    pos[1] = 0;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2 - 1;
                }
            }
            else if (t % 2 == 0) { // cos_sin
                pos[0] = (mode_s + mode_t) * 2;
                if (mode_s == mode_t) {
                    pos[1] = -999;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2;
                }
            }
        }

        return pos;
    };

    std::function<double*(int, int, int)> coefficients_v_ex = [&] (int s, int t, int nfp_ft) {
        double *coeff = new double [2]();
        double sign = -999999;
        double factor;

        if (s == 0) {
            factor = 0.0;
        }
        else if (s % 2 == 1) {
            factor = - (s+1)/2 * nfp_ft;
        }
        else if (s % 2 == 0) {
            factor =   (s+0)/2 * nfp_ft;
        }
        
        if ((s - t) != 0) {
            sign = (s - t) / abs(s - t);
        }

        if (s == 0 || t == 0) { // 0_X X_0
            coeff[0] = 0.5;
            coeff[1] = 0.5;
        }
        else if (s % 2 == 1) {
            if (t % 2 == 1) { // sin_cos
                coeff[0] =  0.5;
                coeff[1] =  0.5 * sign;
            }
            else if (t % 2 == 0) { // sin_sin
                coeff[0] = -0.5;
                coeff[1] =  0.5;
            }
        }
        else if (s % 2 == 0) {
            if (t % 2 == 1) { // cos_cos
                coeff[0] =  0.5;
                coeff[1] =  0.5;
            }
            else if (t % 2 == 0) { // cos_sin
                coeff[0] =  0.5;
                coeff[1] = -0.5 * sign;
            }
        }

        coeff[0] *= factor;
        coeff[1] *= factor;
        
        return coeff;
    };

    vector<int> vec3;
    int flag3 = 0;
    double val1 = 0, val2 = 0;
    // divergence-free constraint
    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {

            for (int k = 0; k < vec3.size(); k++) {
                if (vec3[k] == points_sequence[i][2+3*j]) {
                    flag3 = 1;
                    break;
                }
            }
            if (flag3 == 1) {
                flag3 = 0;
                continue;
            }
            else {
                vec3.push_back(points_sequence[i][2+3*j]);
            }

            index2 = points_sequence[i][2+3*j] * num_modes_divergence;
            index3 = points_sequence[i][2+3*j] * (6*num_modes_scalar);
            for (int s = 0; s < num_modes; s++) { // mg
                for (int t = 0; t < num_modes; t++) { // coeff
                    for (int k = 0; k < num_variables; k++) {
                        if (k != 2) {
                            p1 = position_s(s, t);
                            p2 = coefficients_v(s, t, nfp_ft);
                        }
                        else {
                            p1 = position_s_ex(s,t);
                            p2 = coefficients_v_ex(s, t, nfp_ft);
                        }
                        index1 = (k / num_scalar) * num_modes_scalar + s * num_scalar + (k % num_scalar);
                        for (int l = 0; l < 2; l++) {
                            if (p1[l] == -999 || p1[l] >= num_modes_divergence) {
                                continue;
                            }
                            else if (p1[l] < 0) { //  || p1[l] >= num_modes_divergence
                                printf("k = %d, p1[l] = %d, exit...\n", k, p1[l]);
                                fflush(stdout);
                                exit(0);
                            }
                            // coeff_A.push_back( T(num_extended_DoFs + index2 + p1[l], index3 + index1, divergence_tri_3_var_modes[i][j][k][t] * p2[l]) );
                            // coeff_A.push_back( T(index3 + index1, num_extended_DoFs + index2 + p1[l], divergence_tri_3_var_modes[i][j][k][t] * p2[l]) );
                        }
                        delete [] p1;
                        delete [] p2;
                    }
                }
            }

        }
    }

    A.setFromTriplets(coeff_A.begin(), coeff_A.end());
    A.makeCompressed();
    vector<T>().swap(coeff_A);


    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            index1 = 6*num_modes_scalar*j;
            index2 = points_sequence[i][2+3*j]*6*num_modes_scalar;
            for (int s = 0; s < 6*num_modes_scalar; s++) {
                B(index2 + s) += B_tri_18Nn[i][index1 + s];
            }

            if (points_sequence[i][0+3*j] == 1) {
                index2 = points_sequence[i][1+3*j]*num_modes_scalar;
                index3 = points_sequence[i][1+3*j]*num_modes;
                for (int s = 0; s < num_modes; s++) {
                    for (int t = 0; t < num_scalar; t++) { // num_edges*num_modes
                        // B(num_total_DoFs + index2 + num_scalar*s + t) = boundary_DoFs[t][index3 + s];
                    }
                }
            }

        }
    }

    timing(time_s);

    printf("solving linear equation...\n");
    fflush(stdout);    

    printf(" preconditioning...\n");
    solver2.compute(A);
    if (solver2.info() != Eigen::Success) {
        printf("solve_magnetic_field_sparse_lagrange_multiplier: decomposition failed...\n");
        exit(0);
    }
    fflush(stdout); 

    // solver2.setMaxIterations(1999);
    X = solver2.solve(B);
    // printf("estimated error = %0.9f, iterations = %d\n", solver2.error(), solver2.iterations());
    printf("estimated error = %0.9f\n", abs((A * X - B).norm()) / (B.norm()));
    timing(time_s);

    for (int i = 0; i < num_total_DoFs; i++) {
        X_upper(i) = X(i);
    }


    printf(" renumbering...\n");
    fflush(stdout);
    for (int i = 0; i < num_total_DoFs; i++) {
        X_re[i%num_scalar][i/num_scalar] = X_upper(i);
    }

    for (int j = 0; j < num_scalar; j++) {
        for (int i = 0; i < num_total_DoFs/num_scalar; i++) {
            printf("%0.9f  ", X_re[j][i]);
            if ((i+1)%(num_modes) == 0) {
                printf("        ");
            }
            if ((i+1)%(6*num_modes) == 0) {
                printf("\n");
            }
        }
        printf("\n");
    }
    printf("\n");
    fflush(stdout);
    // exit(0);

    int kk = 0;
    double ***F_sca_tri_18N = new double **[num_scalar];
    for (int i = 0; i < num_scalar; i++) {
        F_sca_tri_18N[i] = new double *[num_triangels_in_plane];
        for (int j = 0; j < num_triangels_in_plane; j++) {
            F_sca_tri_18N[i][j] = new double [18*num_modes] ();
            for (int k = 0; k < 18*num_modes; k++) {
                kk = k/(6*num_modes);
                F_sca_tri_18N[i][j][k] = X_re[i][points_sequence[j][2+3*kk]*6*num_modes + k - kk*6*num_modes];
            }
        }
    } 

    delete [] X_re;

    return F_sca_tri_18N;

}

double ***solve_magnetic_field_sparse_ALM_div(double**** divergence_tri_3_var_modes, double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, double **normal_DoFs,  int num_modes, int num_triangels_in_plane, int nfp_ft, int num_scalar, int flag, double error, time_t time_s) {

    printf("stiffness matrix assembling...\n");
    fflush(stdout);

    int num_variables = 9;
    int num_edges = points_sequence[0][9];
    int num_vertices = points_sequence[0][10];
    int num_modes_divergence = num_modes;
    // int num_modes_divergence = 0;
    int num_modes_normal = num_modes;

    int num_modes_scalar = num_modes*num_scalar;
    int num_total_DoFs = num_vertices * (6*num_modes_scalar);
    int num_divergence_DoFs = num_vertices * num_modes_divergence;
    int num_boundary_DoFs;
    if (flag == 0) {
        num_boundary_DoFs = num_edges * num_modes_scalar;
    }
    else if (flag == 1) {
        num_boundary_DoFs = num_edges * num_modes_normal + 1;
    }
    else {
        printf("dirichlet condition needed...\n");
        exit(0);
    }
    
    int num_extended_DoFs = num_boundary_DoFs + num_divergence_DoFs;
    
    SparseMatrix<double> A(num_total_DoFs, num_total_DoFs);
    A.setZero();

    SparseMatrix<double> L(num_extended_DoFs, num_total_DoFs);
    L.setZero();

    VectorXd B(num_total_DoFs);
    B.setZero();

    VectorXd BL(num_extended_DoFs);
    BL.setZero();

    VectorXd tau(num_extended_DoFs);
    tau.setOnes();

    VectorXd X(num_total_DoFs);
    VectorXd X_upper(num_total_DoFs);
    X.setZero();
    X_upper.setZero();

    double **X_re = new double *[num_scalar];
    for (int i = 0; i < num_scalar; i++) {
        X_re[i] = new double [num_vertices*6*num_modes]();
    }

    typedef Eigen::Triplet<double> T;
    vector<T> coeff_A, coeff_L;

    // Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver2;
    // Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver1;
    // Eigen::LeastSquaresConjugateGradient<SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver1;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver2;
    // Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver2;

    int index_J, index_K, index_j1, index_j2, index_k1;

    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            index_J = 6*num_modes_scalar*j;
            index_j1 = points_sequence[i][2+3*j]*6*num_modes_scalar;
            for (int k = 0; k < 3; k++) {
                index_K = 6*num_modes_scalar*k;
                index_k1 = points_sequence[i][2+3*k]*6*num_modes_scalar;
                for (int s = 0; s < 6*num_modes_scalar; s++) {
                    for (int t = 0; t < 6*num_modes_scalar; t++) {
                        coeff_A.push_back( T(index_j1 + s, index_k1 + t, M_tri_18Nn_18Nn[i][index_J + s][index_K + t]) );
                    }
                }
            }
        }
    }

    
    std::function<int*(int, int)> position_s = [&] (int s, int t) {
        int *pos = new int [2]();

        int mode_s;
        int mode_t;
        if (s % 2 == 1) {
            mode_s = (s + 1) / 2;
        }
        else {
            mode_s = s / 2;
        }
        if (t % 2 == 1) {
            mode_t = (t + 1) / 2;
        }
        else {
            mode_t = t / 2;
        }

        if (s == 0) { // 0_X
            pos[0] = t;
            pos[1] = t;
        }
        else if (t == 0) { // X_0
            pos[0] = s;
            pos[1] = s;
        }
        else if (s % 2 == 1) {
            if (t % 2 == 1) { // cos_cos
                pos[0] = (mode_s + mode_t) * 2 - 1;
                if (mode_s == mode_t) {
                    pos[1] = 0;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2 - 1;
                }
            }
            else if (t % 2 == 0) { // cos_sin
                pos[0] = (mode_s + mode_t) * 2;
                if (mode_s == mode_t) {
                    pos[1] = -999;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2;
                }
            }
        }
        else if (s % 2 == 0) {
            if (t % 2 == 1) { // sin_cos
                pos[0] = (mode_s + mode_t) * 2;
                if (mode_s == mode_t) {
                    pos[1] = -999;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2;
                }
            }
            else if (t % 2 == 0) { // sin_sin
                pos[0] = (mode_s + mode_t) * 2 - 1;
                if (mode_s == mode_t) {
                    pos[1] = 0;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2 - 1;
                }
            }
        }

        return pos;
    };

    std::function<double*(int, int, int)> coefficients_v = [&] (int s, int t, int nfp_ft) {
        double *coeff = new double [2]();

        double sign = -999999;
        if ((s - t) != 0) {
            sign = (s - t) / abs(s - t);
        }

        if (s == 0 || t == 0) { // 0_X X_0
            coeff[0] = 0.5;
            coeff[1] = 0.5;
        }
        else if (s % 2 == 1) {
            if (t % 2 == 1) { // cos_cos
                coeff[0] =  0.5;
                coeff[1] =  0.5;
            }
            else if (t % 2 == 0) { // cos_sin
                coeff[0] =  0.5;
                coeff[1] = -0.5 * sign;
            }
        }
        else if (s % 2 == 0) {
            if (t % 2 == 1) { // sin_cos
                coeff[0] =  0.5;
                coeff[1] =  0.5 * sign;
            }
            else if (t % 2 == 0) { // sin_sin
                coeff[0] = -0.5;
                coeff[1] =  0.5;
            }
        }

        return coeff;
    };

    int *p1;
    double *p2;
    vector<int> vec2;
    int flag2 = 0;
    int index1, index2, index3;
    int boundary_point;
    // dirichlet boundary condition 
    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {    
            if (points_sequence[i][0+3*j] == 1) {

                if (points_sequence[i][1+3*j] == 0) {
                    boundary_point = points_sequence[i][2+3*j];
                }

                for (int k = 0; k < vec2.size(); k++) {
                    if (vec2[k] == points_sequence[i][1+3*j]) {
                        flag2 = 1;
                        break;
                    }
                }
                if (flag2 == 1) {
                    flag2 = 0;
                    continue;
                }
                else {
                    vec2.push_back(points_sequence[i][1+3*j]);
                }

                if (flag == 0) {
                    index2 = points_sequence[i][1+3*j]*num_modes_scalar;
                    index3 = points_sequence[i][2+3*j]*6*num_modes_scalar;
                    for (int s = 0; s < num_modes; s++) {
                        for (int t = 0; t < num_scalar; t++) {
                            coeff_L.push_back( T(index2 + num_scalar*s + t, index3 + num_scalar*s + t, 1.0) );
                        }
                    }
                }
                else if (flag == 1) {
                    index2 = points_sequence[i][1+3*j] * num_modes_normal;
                    index3 = points_sequence[i][2+3*j] * (6*num_modes_scalar);
                    index1 = points_sequence[i][1+3*j] * (2*num_modes-1);
                    for (int s = 0; s < num_modes; s++) { // mag
                        for (int t = 0; t < 2*num_modes-1; t++) { // coeff
                            p1 = position_s(s, t);
                            p2 = coefficients_v(s, t, nfp_ft);
                            for (int k = 0; k < 3; k++) {
                                for (int l = 0; l < 2; l++) {
                                    if (p1[l] == -999 || p1[l] >= num_modes_normal) {
                                        continue;
                                    }
                                    coeff_L.push_back( T(index2 + p1[l], index3 + s * num_scalar + k, normal_DoFs[k][index1 + t] * p2[l]) );
                                }
                            }
                            delete [] p1;
                            delete [] p2;
                        }
                    }
                }
                
            }

        }
    }

    if (flag == 1) {
        // BL(num_boundary_DoFs - 2) = 1.0;
        // coeff_L.push_back( T(num_boundary_DoFs - 2, boundary_point * (6*num_modes_scalar) + 2, 1) );
        // BL(num_boundary_DoFs - 1) = 0.1;
        // coeff_L.push_back( T(num_boundary_DoFs - 1, boundary_point * (6*num_modes_scalar) + 1, 1) );

        BL(num_boundary_DoFs - 1) = 0.85;
        coeff_L.push_back( T(num_boundary_DoFs - 1, boundary_point * (6*num_modes_scalar) + 2, 1) );
    }



    std::function<int*(int, int)> position_s_ex = [&] (int s, int t) {
        int *pos = new int [2]();

        int mode_s;
        int mode_t;
        if (s % 2 == 1) {
            mode_s = (s + 1) / 2;
        }
        else {
            mode_s = s / 2;
        }
        if (t % 2 == 1) {
            mode_t = (t + 1) / 2;
        }
        else {
            mode_t = t / 2;
        }
        
        if (s == 0) { // 0_X
            pos[0] = -999;
            pos[1] = -999;
        }
        else if (t == 0) { // X_0
            if (s % 2 == 1) {
                pos[0] = s + 1;
                pos[1] = s + 1;
            }
            else if (s % 2 == 0) {
                pos[0] = s - 1;
                pos[1] = s - 1;
            }
        }
        else if (s % 2 == 1) {
            if (t % 2 == 1) { // sin_cos
                pos[0] = (mode_s + mode_t) * 2;
                if (mode_s == mode_t) {
                    pos[1] = -999;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2;
                }
            }
            else if (t % 2 == 0) { // sin_sin
                pos[0] = (mode_s + mode_t) * 2 - 1;
                if (mode_s == mode_t) {
                    pos[1] = 0;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2 - 1;
                }
            }
        }
        else if (s % 2 == 0) {
            if (t % 2 == 1) { // cos_cos
                pos[0] = (mode_s + mode_t) * 2 - 1;
                if (mode_s == mode_t) {
                    pos[1] = 0;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2 - 1;
                }
            }
            else if (t % 2 == 0) { // cos_sin
                pos[0] = (mode_s + mode_t) * 2;
                if (mode_s == mode_t) {
                    pos[1] = -999;
                }
                else {
                    pos[1] = abs(mode_s - mode_t) * 2;
                }
            }
        }

        return pos;
    };

    std::function<double*(int, int, int)> coefficients_v_ex = [&] (int s, int t, int nfp_ft) {
        double *coeff = new double [2]();
        double sign = -999999;
        double factor;

        if (s == 0) {
            factor = 0.0;
        }
        else if (s % 2 == 1) {
            factor = - (s+1)/2 * nfp_ft;
        }
        else if (s % 2 == 0) {
            factor =   (s+0)/2 * nfp_ft;
        }
        
        if ((s - t) != 0) {
            sign = (s - t) / abs(s - t);
        }

        if (s == 0 || t == 0) { // 0_X X_0
            coeff[0] = 0.5;
            coeff[1] = 0.5;
        }
        else if (s % 2 == 1) {
            if (t % 2 == 1) { // sin_cos
                coeff[0] =  0.5;
                coeff[1] =  0.5 * sign;
            }
            else if (t % 2 == 0) { // sin_sin
                coeff[0] = -0.5;
                coeff[1] =  0.5;
            }
        }
        else if (s % 2 == 0) {
            if (t % 2 == 1) { // cos_cos
                coeff[0] =  0.5;
                coeff[1] =  0.5;
            }
            else if (t % 2 == 0) { // cos_sin
                coeff[0] =  0.5;
                coeff[1] = -0.5 * sign;
            }
        }

        coeff[0] *= factor;
        coeff[1] *= factor;
        
        return coeff;
    };

    vector<int> vec3;
    int flag3 = 0;
    // divergence-free constraint
    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < vec3.size(); k++) {
                if (vec3[k] == points_sequence[i][2+3*j]) {
                    flag3 = 1;
                    break;
                }
            }
            if (flag3 == 1) {
                flag3 = 0;
                continue;
            }
            else {
                vec3.push_back(points_sequence[i][2+3*j]);
            }
            index2 = points_sequence[i][2+3*j] * num_modes_divergence;
            index3 = points_sequence[i][2+3*j] * (6*num_modes_scalar);
            for (int s = 0; s < num_modes; s++) { // mg
                for (int t = 0; t < 2*num_modes-1; t++) { // coeff
                    for (int k = 0; k < num_variables; k++) {
                        if (k != 2) {
                            p1 = position_s(s, t);
                            p2 = coefficients_v(s, t, nfp_ft);
                        }
                        else {
                            p1 = position_s_ex(s,t);
                            p2 = coefficients_v_ex(s, t, nfp_ft);
                        }
                        index1 = (k / num_scalar) * num_modes_scalar + s * num_scalar + (k % num_scalar);
                        for (int l = 0; l < 2; l++) {
                            if (p1[l] == -999 || p1[l] >= num_modes_divergence) {
                                continue;
                            }
                            else if (p1[l] < 0) { //  || p1[l] >= num_modes_divergence
                                printf("k = %d, p1[l] = %d, exit...\n", k, p1[l]);
                                fflush(stdout);
                                exit(0);
                            }
                            coeff_L.push_back( T(num_boundary_DoFs + index2 + p1[l], index3 + index1, divergence_tri_3_var_modes[i][j][k][t] * p2[l]) );
                        }
                        delete [] p1;
                        delete [] p2;
                    }
                }
            }
        }
    }


    A.setFromTriplets(coeff_A.begin(), coeff_A.end());
    A.makeCompressed();
    vector<T>().swap(coeff_A);

    L.setFromTriplets(coeff_L.begin(), coeff_L.end());
    L.makeCompressed();
    vector<T>().swap(coeff_L);


    for (int i = 0; i < num_triangels_in_plane; i++) {
        for (int j = 0; j < 3; j++) {
            index1 = 6*num_modes_scalar*j;
            index2 = points_sequence[i][2+3*j]*6*num_modes_scalar;
            for (int s = 0; s < 6*num_modes_scalar; s++) {
                B(index2 + s) += B_tri_18Nn[i][index1 + s];
            }

            if (points_sequence[i][0+3*j] == 1 && flag == 0) {
                index2 = points_sequence[i][1+3*j]*num_modes_scalar;
                index3 = points_sequence[i][1+3*j]*num_modes;
                for (int s = 0; s < num_modes; s++) {
                    for (int t = 0; t < num_scalar; t++) { // num_edges*num_modes
                        BL(index2 + num_scalar*s + t) = boundary_DoFs[t][index3 + s];
                    }
                }
            }

        }
    }

    timing(time_s);

    printf("solving linear equation...\n");
    fflush(stdout);    

    printf(" preconditioning...\n");
    double lambda = 1E+03; // 1E+03;
    double rho = 10; // 10;

    for (int i = 0; i < 10; i++) {

        solver2.compute( (A + lambda * (L.transpose() * L)) );
        if (solver2.info() != Eigen::Success) {
            printf("solve_magnetic_field_sparse_lagrange_multiplier: decomposition failed...\n");
            fflush(stdout); 
            exit(0);
        }
        
        X = solver2.solve( (B + lambda * (L.transpose() * BL) - (L.transpose() * tau)) );

        if (lambda < 0.9999999E+07 ) {
            lambda *= rho;
        }
        tau = tau.eval() + lambda * (L * X - BL);

        printf("|L * X - BL| / |BL| = %0.12f\n", (L * X - BL).norm() / (BL).norm());
        fflush(stdout);

    }

    // solver2.compute(A);
    // if (solver2.info() != Eigen::Success) {
    //     printf("solve_magnetic_field_sparse_lagrange_multiplier: decomposition failed...\n");
    //     exit(0);
    // }
    // fflush(stdout); 

    // // solver2.setMaxIterations(1999);
    // X = solver2.solve(B);
    // // printf("estimated error = %0.9f, iterations = %d\n", solver2.error(), solver2.iterations());
    // printf("estimated error = %0.9f\n", abs((A * X - B).norm()) / (B.norm()));
    // timing(time_s);

    for (int i = 0; i < num_total_DoFs; i++) {
        X_upper(i) = X(i);
    }


    printf(" renumbering...\n");
    fflush(stdout);
    for (int i = 0; i < num_total_DoFs; i++) {
        X_re[i%num_scalar][i/num_scalar] = X_upper(i);
    }

    // for (int j = 0; j < num_scalar; j++) {
    //     for (int i = 0; i < num_total_DoFs/num_scalar; i++) {
    //         printf("%0.9f  ", X_re[j][i]);
    //         if ((i+1)%(num_modes) == 0) {
    //             printf("        ");
    //         }
    //         if ((i+1)%(6*num_modes) == 0) {
    //             printf("\n");
    //         }
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    // fflush(stdout);
    // // exit(0);

    int kk = 0;
    double ***F_sca_tri_18N = new double **[num_scalar];
    for (int i = 0; i < num_scalar; i++) {
        F_sca_tri_18N[i] = new double *[num_triangels_in_plane];
        for (int j = 0; j < num_triangels_in_plane; j++) {
            F_sca_tri_18N[i][j] = new double [18*num_modes] ();
            for (int k = 0; k < 18*num_modes; k++) {
                kk = k/(6*num_modes);
                F_sca_tri_18N[i][j][k] = X_re[i][points_sequence[j][2+3*kk]*6*num_modes + k - kk*6*num_modes];
            }
        }
    } 

    delete [] X_re;

    return F_sca_tri_18N;

}









