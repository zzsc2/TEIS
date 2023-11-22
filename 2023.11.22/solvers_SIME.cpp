#include "solvers_SIME.hpp"
#include "superstructure.hpp"

double *FuncJperpendicular(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***B_tri_3_18, double ***gradP_tri_3_18, int position_sec, int position_tri) {

    int dim = 3;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H;
    double uv[20], uv_u[20], uv_v[20];
    double R_ = 0, R_u = 0, R_v = 0, R_zeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_zeta = 0;

    double BR = 0, BZ = 0, Bphi = 0;
    double gradPR = 0, gradPZ = 0, gradPphi = 0;
    double B2;

    for (int i = 0; i < 20; i++) {
        uv[i] =                 power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] * power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] * power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);
        for (int j = 0; j < 18; j++) {
            R_ +=       G[position_tri][i][j] * uv[i] *   R[position_sec][position_tri][j];
            R_u +=      G[position_tri][i][j] * uv_u[i] * R[position_sec][position_tri][j];
            R_v +=      G[position_tri][i][j] * uv_v[i] * R[position_sec][position_tri][j];
            
            Z_ +=       G[position_tri][i][j] * uv[i] *   Z[position_sec][position_tri][j];
            Z_u +=      G[position_tri][i][j] * uv_u[i] * Z[position_sec][position_tri][j];
            Z_v +=      G[position_tri][i][j] * uv_v[i] * Z[position_sec][position_tri][j];

            R_zeta +=   G[position_tri][i][j] * uv[i] *   Rzeta[position_sec][position_tri][j];

            Z_zeta +=   G[position_tri][i][j] * uv[i] *   Zzeta[position_sec][position_tri][j];

            BR +=       G[position_tri][i][j] * uv[i] *   B_tri_3_18[position_tri][0][j];
            BZ +=       G[position_tri][i][j] * uv[i] *   B_tri_3_18[position_tri][1][j];
            Bphi +=     G[position_tri][i][j] * uv[i] *   B_tri_3_18[position_tri][2][j];

            gradPR +=   G[position_tri][i][j] * uv[i] *   gradP_tri_3_18[position_tri][0][j];
            gradPZ +=   G[position_tri][i][j] * uv[i] *   gradP_tri_3_18[position_tri][1][j];
            gradPphi += G[position_tri][i][j] * uv[i] *   gradP_tri_3_18[position_tri][2][j];
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

    B2 = BR * BR + BZ * BZ + Bphi * Bphi;

    // R  Z  φ
    // 0  1  2

    double *var1 = new double [3]();

    var1[0] = (BZ * gradPphi - Bphi * gradPZ) / B2;
    var1[1] = (Bphi * gradPR - BR * gradPphi) / B2;
    var1[2] = (BR * gradPZ   - BZ * gradPR) / B2;

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}

double *FuncGradSigma(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****Jperpendicular_2_tri_3_18, double ***B_tri_3_18, int position_sec, int position_tri) {

    int dim = 3;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H;
    double uv[20], uv_u[20], uv_v[20];
    double R_ = 0, R_u = 0, R_v = 0, R_zeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_zeta = 0;

    double BR = 0, BZ = 0, Bphi = 0;
    double J[3] = {0};
    double J_u[3] = {0}, J_v[3] = {0}, J_zeta[3] = {0};
    double J_R[3], J_Z[3], J_phi[3];

    for (int i = 0; i < 20; i++) {
        uv[i] =                 power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] * power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] * power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);
        for (int j = 0; j < 18; j++) {
            R_ +=            G[position_tri][i][j] * uv[i] *   R[position_sec][position_tri][j];
            R_u +=           G[position_tri][i][j] * uv_u[i] * R[position_sec][position_tri][j];
            R_v +=           G[position_tri][i][j] * uv_v[i] * R[position_sec][position_tri][j];
            
            Z_ +=            G[position_tri][i][j] * uv[i] *   Z[position_sec][position_tri][j];
            Z_u +=           G[position_tri][i][j] * uv_u[i] * Z[position_sec][position_tri][j];
            Z_v +=           G[position_tri][i][j] * uv_v[i] * Z[position_sec][position_tri][j];

            R_zeta +=        G[position_tri][i][j] * uv[i] *   Rzeta[position_sec][position_tri][j];

            Z_zeta +=        G[position_tri][i][j] * uv[i] *   Zzeta[position_sec][position_tri][j];

            BR +=            G[position_tri][i][j] * uv[i] *   B_tri_3_18[position_tri][0][j];
            BZ +=            G[position_tri][i][j] * uv[i] *   B_tri_3_18[position_tri][1][j];
            Bphi +=          G[position_tri][i][j] * uv[i] *   B_tri_3_18[position_tri][2][j];

            for (int k = 0; k < 3; k++) {
                J[k] +=      G[position_tri][i][j] * uv[i] *   Jperpendicular_2_tri_3_18[0][position_tri][k][j];
                J_u[k] +=    G[position_tri][i][j] * uv_u[i] * Jperpendicular_2_tri_3_18[0][position_tri][k][j];
                J_v[k] +=    G[position_tri][i][j] * uv_v[i] * Jperpendicular_2_tri_3_18[0][position_tri][k][j];
                J_zeta[k] += G[position_tri][i][j] * uv[i] *   Jperpendicular_2_tri_3_18[1][position_tri][k][j];
            }
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

    
    for (int i = 0; i < 3; i++) {
        J_R[i] =   p[0][0] * J_u[i] + p[0][1] * J_v[i] + p[0][2] * J_zeta[i];
        J_Z[i] =   p[1][0] * J_u[i] + p[1][1] * J_v[i] + p[1][2] * J_zeta[i];
        J_phi[i] = p[2][0] * J_u[i] + p[2][1] * J_v[i] + p[2][2] * J_zeta[i];
    }

    double B2 = BR * BR + BZ * BZ + Bphi * Bphi;
    double divJperpendicular = (J[0] / R_) + (J_R[0]) + (J_Z[1]) + (J_phi[2] / R_);

    // R  Z  φ
    // 0  1  2

    double *var1 = new double [3]();

    var1[0] = - divJperpendicular / B2 * BR;
    var1[1] = - divJperpendicular / B2 * BZ;
    var1[2] = - divJperpendicular / B2 * Bphi; 

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}

double FuncDivGradSigma(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****GradSigma_2_tri_3_18, int position_sec, int position_tri) {
    
    int dim = 3;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H;
    double uv[20], uv_u[20], uv_v[20];
    double R_ = 0, R_u = 0, R_v = 0, R_zeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_zeta = 0;

    double GradSigma[3] = {0};
    double GradSigma_u[3] = {0}, GradSigma_v[3] = {0}, GradSigma_zeta[3] = {0};
    double GradSigma_R[3], GradSigma_Z[3], GradSigma_phi[3];

    for (int i = 0; i < 20; i++) {
        uv[i] =                 power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] * power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] * power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);
        for (int j = 0; j < 18; j++) {
            R_ +=            G[position_tri][i][j] * uv[i] *   R[position_sec][position_tri][j];
            R_u +=           G[position_tri][i][j] * uv_u[i] * R[position_sec][position_tri][j];
            R_v +=           G[position_tri][i][j] * uv_v[i] * R[position_sec][position_tri][j];
            
            Z_ +=            G[position_tri][i][j] * uv[i] *   Z[position_sec][position_tri][j];
            Z_u +=           G[position_tri][i][j] * uv_u[i] * Z[position_sec][position_tri][j];
            Z_v +=           G[position_tri][i][j] * uv_v[i] * Z[position_sec][position_tri][j];

            R_zeta +=        G[position_tri][i][j] * uv[i] *   Rzeta[position_sec][position_tri][j];

            Z_zeta +=        G[position_tri][i][j] * uv[i] *   Zzeta[position_sec][position_tri][j];

            for (int k = 0; k < 3; k++) {
                GradSigma[k] +=      G[position_tri][i][j] * uv[i] *   GradSigma_2_tri_3_18[0][position_tri][k][j];
                GradSigma_u[k] +=    G[position_tri][i][j] * uv_u[i] * GradSigma_2_tri_3_18[0][position_tri][k][j];
                GradSigma_v[k] +=    G[position_tri][i][j] * uv_v[i] * GradSigma_2_tri_3_18[0][position_tri][k][j];
                GradSigma_zeta[k] += G[position_tri][i][j] * uv[i] *   GradSigma_2_tri_3_18[1][position_tri][k][j];
            }
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
    
    // R  Z  φ
    // 0  1  2

    for (int i = 0; i < 3; i++) {
        GradSigma_R[i] =   p[0][0] * GradSigma_u[i] + p[0][1] * GradSigma_v[i] + p[0][2] * GradSigma_zeta[i];
        GradSigma_Z[i] =   p[1][0] * GradSigma_u[i] + p[1][1] * GradSigma_v[i] + p[1][2] * GradSigma_zeta[i];
        GradSigma_phi[i] = p[2][0] * GradSigma_u[i] + p[2][1] * GradSigma_v[i] + p[2][2] * GradSigma_zeta[i];
    }

    double var1 = ((GradSigma[0] / R_) + (GradSigma_R[0]) + (GradSigma_Z[1]) + (GradSigma_phi[2] / R_)) * D * R_;

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}

double *FuncJparallel(double u, double v, double ***G, double ***B_tri_3_18, double **sigma_tri_18, int position_sec, int position_tri) {

    double uv[20];

    double BR = 0, BZ = 0, Bphi = 0, sigma = 0;
    double pressure_u = 0, pressure_v = 0, pressure_zeta = 0;
    double B2, pressure_R, pressure_Z, pressure_phi;

    for (int i = 0; i < 20; i++) {
        uv[i] = power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        for (int j = 0; j < 18; j++) {
            BR +=     G[position_tri][i][j] * uv[i] * B_tri_3_18[position_tri][0][j];
            BZ +=     G[position_tri][i][j] * uv[i] * B_tri_3_18[position_tri][1][j];
            Bphi +=   G[position_tri][i][j] * uv[i] * B_tri_3_18[position_tri][2][j];

            sigma +=  G[position_tri][i][j] * uv[i] * sigma_tri_18[position_tri][j];
        }
    }

    double *var1 = new double [3]();

    var1[0] = sigma * BR;
    var1[1] = sigma * BZ;
    var1[2] = sigma * Bphi;

    return var1;
}

double *FuncCurlJ(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****J_2_tri_3_18, int position_sec, int position_tri) {

    int dim = 3;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H;
    double uv[20], uv_u[20], uv_v[20];
    double R_ = 0, R_u = 0, R_v = 0, R_zeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_zeta = 0;

    double J[3] = {0,0,0};
    double J_u[3] = {0,0,0}, J_v[3] = {0,0,0}, J_zeta[3] = {0,0,0};
    double J_R[3], J_Z[3], J_phi[3];

    for (int i = 0; i < 20; i++) {
        uv[i] =                 power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] * power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] * power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);
        for (int j = 0; j < 18; j++) {
            R_ +=            G[position_tri][i][j] * uv[i] *   R[position_sec][position_tri][j];
            R_u +=           G[position_tri][i][j] * uv_u[i] * R[position_sec][position_tri][j];
            R_v +=           G[position_tri][i][j] * uv_v[i] * R[position_sec][position_tri][j];
            
            Z_ +=            G[position_tri][i][j] * uv[i] *   Z[position_sec][position_tri][j];
            Z_u +=           G[position_tri][i][j] * uv_u[i] * Z[position_sec][position_tri][j];
            Z_v +=           G[position_tri][i][j] * uv_v[i] * Z[position_sec][position_tri][j];

            R_zeta +=        G[position_tri][i][j] * uv[i] *   Rzeta[position_sec][position_tri][j];

            Z_zeta +=        G[position_tri][i][j] * uv[i] *   Zzeta[position_sec][position_tri][j];

            for (int k = 0; k < 3; k++) {
                J[k] +=      G[position_tri][i][j] * uv[i] *   J_2_tri_3_18[0][position_tri][k][j];
                J_u[k] +=    G[position_tri][i][j] * uv_u[i] * J_2_tri_3_18[0][position_tri][k][j];
                J_v[k] +=    G[position_tri][i][j] * uv_v[i] * J_2_tri_3_18[0][position_tri][k][j];
                J_zeta[k] += G[position_tri][i][j] * uv[i] *   J_2_tri_3_18[1][position_tri][k][j];
            }
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

    for (int i = 0; i < 3; i++) {
        J_R[i] =   p[0][0] * J_u[i] + p[0][1] * J_v[i] + p[0][2] * J_zeta[i];
        J_Z[i] =   p[1][0] * J_u[i] + p[1][1] * J_v[i] + p[1][2] * J_zeta[i];
        J_phi[i] = p[2][0] * J_u[i] + p[2][1] * J_v[i] + p[2][2] * J_zeta[i];
    }

    // R  Z  φ
    // 0  1  2

    double *var1 = new double [3]();

    var1[0] = - mu0 * (J_phi[1] / R_ - J_Z[2]) * D * R_;     // JZ_phi / R - Jphi_Z
    var1[1] = - mu0 * (J[2] / R_ + J_R[2] - J_phi[0] / R_) * D * R_; // Jphi / R + Jphi_R - JR_phi / R
    var1[2] = - mu0 * (J_Z[0] - J_R[1]) * D * R_;    // JR_Z - JZ_R

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}



double *FuncGradP(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***pressure_2_tri_18, int position_sec, int position_tri) {

    int dim = 3;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H;
    double uv[20], uv_u[20], uv_v[20];
    double R_ = 0, R_u = 0, R_v = 0, R_zeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_zeta = 0;

    double BR = 0, BZ = 0, Bphi = 0;
    double pressure_u = 0, pressure_v = 0, pressure_zeta = 0;
    double B2, pressure_R, pressure_Z, pressure_phi;

    for (int i = 0; i < 20; i++) {
        uv[i] =                 power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] * power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] * power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);
        for (int j = 0; j < 18; j++) {
            R_ +=            G[position_tri][i][j] * uv[i] *   R[position_sec][position_tri][j];
            R_u +=           G[position_tri][i][j] * uv_u[i] * R[position_sec][position_tri][j];
            R_v +=           G[position_tri][i][j] * uv_v[i] * R[position_sec][position_tri][j];
            
            Z_ +=            G[position_tri][i][j] * uv[i] *   Z[position_sec][position_tri][j];
            Z_u +=           G[position_tri][i][j] * uv_u[i] * Z[position_sec][position_tri][j];
            Z_v +=           G[position_tri][i][j] * uv_v[i] * Z[position_sec][position_tri][j];

            R_zeta +=        G[position_tri][i][j] * uv[i] *   Rzeta[position_sec][position_tri][j];

            Z_zeta +=        G[position_tri][i][j] * uv[i] *   Zzeta[position_sec][position_tri][j];

            pressure_u +=    G[position_tri][i][j] * uv_u[i] * pressure_2_tri_18[0][position_tri][j];
            pressure_v +=    G[position_tri][i][j] * uv_v[i] * pressure_2_tri_18[0][position_tri][j];
            pressure_zeta += G[position_tri][i][j] * uv[i] *   pressure_2_tri_18[1][position_tri][j];
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

    pressure_R =   p[0][0] * pressure_u + p[0][1] * pressure_v + p[0][2] * pressure_zeta;
    pressure_Z =   p[1][0] * pressure_u + p[1][1] * pressure_v + p[1][2] * pressure_zeta;
    pressure_phi = p[2][0] * pressure_u + p[2][1] * pressure_v + p[2][2] * pressure_zeta;

    // R  Z  φ
    // 0  1  2

    double *var1 = new double [3]();

    var1[0] = pressure_R;
    var1[1] = pressure_Z;
    var1[2] = pressure_phi / R_;

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}

double *FuncGradPparallel(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ***B_tri_3_18, double ***gradP_tri_3_18, int position_sec, int position_tri) {

    int dim = 3;

    double uv[20], uv_u[20], uv_v[20];

    double BR = 0, BZ = 0, Bphi = 0;
    double gradPR = 0, gradPZ = 0, gradPphi = 0;

    for (int i = 0; i < 20; i++) {
        uv[i] =                 power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] * power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] * power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);
        for (int j = 0; j < 18; j++) {
            BR +=       G[position_tri][i][j] * uv[i] *  B_tri_3_18[position_tri][0][j];
            BZ +=       G[position_tri][i][j] * uv[i] *  B_tri_3_18[position_tri][1][j];
            Bphi +=     G[position_tri][i][j] * uv[i] *  B_tri_3_18[position_tri][2][j];

            gradPR +=   G[position_tri][i][j] * uv[i] *  gradP_tri_3_18[position_tri][0][j];
            gradPZ +=   G[position_tri][i][j] * uv[i] *  gradP_tri_3_18[position_tri][1][j];
            gradPphi += G[position_tri][i][j] * uv[i] *  gradP_tri_3_18[position_tri][2][j];
        }
    }

    // R  Z  φ
    // 0  1  2

    double projection = 1.0 * (gradPR * BR + gradPZ * BZ + gradPphi * Bphi) / (BR * BR + BZ * BZ + Bphi * Bphi);

    double *var1 = new double [3]();

    var1[0] = gradPR - projection * BR;
    var1[1] = gradPZ - projection * BZ;
    var1[2] = gradPphi - projection * Bphi;

    // var1[0] = gradPR - 0.999 * projection * BR;
    // var1[1] = gradPZ - 0.999 * projection * BZ;
    // var1[2] = gradPphi - 0.999 * projection * Bphi;

    // var1[0] = gradPR - 0;
    // var1[1] = gradPZ - 0;
    // var1[2] = gradPphi - 0;

    return var1;
}

double FuncDivGradPparallel(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***G, double ****GradPparallel_2_tri_3_18, int position_sec, int position_tri) {

    int dim = 3;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H;
    double uv[20], uv_u[20], uv_v[20];
    double R_ = 0, R_u = 0, R_v = 0, R_zeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_zeta = 0;

    double GradPparallel[3] = {0};
    double GradPparallel_u[3] = {0}, GradPparallel_v[3] = {0}, GradPparallel_zeta[3] = {0};
    double GradPparallel_R[3], GradPparallel_Z[3], GradPparallel_phi[3];

    for (int i = 0; i < 20; i++) {
        uv[i] =                 power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]);
        uv_u[i] =  m_order[i] * power_nonnegative(u, m_order[i]-1) * power_nonnegative(v, n_order[i]);
        uv_v[i] =  n_order[i] * power_nonnegative(u, m_order[i]) *   power_nonnegative(v, n_order[i]-1);
        for (int j = 0; j < 18; j++) {
            R_ +=            G[position_tri][i][j] * uv[i] *   R[position_sec][position_tri][j];
            R_u +=           G[position_tri][i][j] * uv_u[i] * R[position_sec][position_tri][j];
            R_v +=           G[position_tri][i][j] * uv_v[i] * R[position_sec][position_tri][j];
            
            Z_ +=            G[position_tri][i][j] * uv[i] *   Z[position_sec][position_tri][j];
            Z_u +=           G[position_tri][i][j] * uv_u[i] * Z[position_sec][position_tri][j];
            Z_v +=           G[position_tri][i][j] * uv_v[i] * Z[position_sec][position_tri][j];

            R_zeta +=        G[position_tri][i][j] * uv[i] *   Rzeta[position_sec][position_tri][j];

            Z_zeta +=        G[position_tri][i][j] * uv[i] *   Zzeta[position_sec][position_tri][j];

            for (int k = 0; k < 3; k++) {
                GradPparallel[k] +=      G[position_tri][i][j] * uv[i] *   GradPparallel_2_tri_3_18[0][position_tri][k][j];
                GradPparallel_u[k] +=    G[position_tri][i][j] * uv_u[i] * GradPparallel_2_tri_3_18[0][position_tri][k][j];
                GradPparallel_v[k] +=    G[position_tri][i][j] * uv_v[i] * GradPparallel_2_tri_3_18[0][position_tri][k][j];
                GradPparallel_zeta[k] += G[position_tri][i][j] * uv[i] *   GradPparallel_2_tri_3_18[1][position_tri][k][j];
            }
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
    
    // R  Z  φ
    // 0  1  2

    for (int i = 0; i < 3; i++) {
        GradPparallel_R[i] =   p[0][0] * GradPparallel_u[i] + p[0][1] * GradPparallel_v[i] + p[0][2] * GradPparallel_zeta[i];
        GradPparallel_Z[i] =   p[1][0] * GradPparallel_u[i] + p[1][1] * GradPparallel_v[i] + p[1][2] * GradPparallel_zeta[i];
        GradPparallel_phi[i] = p[2][0] * GradPparallel_u[i] + p[2][1] * GradPparallel_v[i] + p[2][2] * GradPparallel_zeta[i];
    }

    double var1 = ((GradPparallel[0] / R_) + (GradPparallel_R[0]) + (GradPparallel_Z[1]) + (GradPparallel_phi[2] / R_)) * D * R_;

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}










double FuncPressure(double u, double v, double ***R, double ***Z, double ***Rzeta, double ***Zzeta, double ***Rzetazeta, double ***Zzetazeta, double ***G, double ***pressure_3_tri_18, double ****GradPparallel_2_tri_3_18, double ***B_tri_3_18, int position_sec, int position_tri) {

    int dim = 9;

    double **p = new double *[dim];
    for (int i = 0; i < dim; i++) {
        p[i] = new double [dim]();
    }

    double D, E, H, alpha[5], beta[5], gamma[5];
    double A_31, A_32, A_41, A_42, A_51, A_52;
    double uv[20], uv_u[20], uv_v[20], uv_uu[20], uv_uv[20], uv_vv[20];
    double R_ = 0, R_u = 0, R_v = 0, R_uu = 0, R_uv = 0, R_vv = 0, R_zetau = 0, R_zetav = 0, R_zetazeta = 0, R_zeta = 0;
    double pres_u = 0, pres_v = 0, pres_uu = 0, pres_uv = 0, pres_vv = 0, pres_zetau = 0, pres_zetav = 0, pres_zetazeta = 0, pres_zeta = 0;
    double Z_ = 0, Z_u = 0, Z_v = 0, Z_uu = 0, Z_uv = 0, Z_vv = 0, Z_zetau = 0, Z_zetav = 0, Z_zetazeta = 0, Z_zeta = 0;

    double pres_R, pres_RR, pres_ZZ, pres_phiphi;

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

            pres_u +=     G[position_tri][i][j] * uv_u[i] *  pressure_3_tri_18[0][position_tri][j];
            pres_v +=     G[position_tri][i][j] * uv_v[i] *  pressure_3_tri_18[0][position_tri][j];
            pres_uu +=    G[position_tri][i][j] * uv_uu[i] * pressure_3_tri_18[0][position_tri][j];
            pres_uv +=    G[position_tri][i][j] * uv_uv[i] * pressure_3_tri_18[0][position_tri][j];
            pres_vv +=    G[position_tri][i][j] * uv_vv[i] * pressure_3_tri_18[0][position_tri][j];

            R_zeta +=     G[position_tri][i][j] * uv[i] *    Rzeta[position_sec][position_tri][j];
            R_zetau +=    G[position_tri][i][j] * uv_u[i] *  Rzeta[position_sec][position_tri][j];
            R_zetav +=    G[position_tri][i][j] * uv_v[i] *  Rzeta[position_sec][position_tri][j];

            Z_zeta +=     G[position_tri][i][j] * uv[i] *    Zzeta[position_sec][position_tri][j];
            Z_zetau +=    G[position_tri][i][j] * uv_u[i] *  Zzeta[position_sec][position_tri][j];
            Z_zetav +=    G[position_tri][i][j] * uv_v[i] *  Zzeta[position_sec][position_tri][j];

            pres_zeta +=  G[position_tri][i][j] * uv[i] *    pressure_3_tri_18[1][position_tri][j];
            pres_zetau += G[position_tri][i][j] * uv_u[i] *  pressure_3_tri_18[1][position_tri][j];
            pres_zetav += G[position_tri][i][j] * uv_v[i] *  pressure_3_tri_18[1][position_tri][j];

            R_zetazeta += G[position_tri][i][j] * uv[i] *    Rzetazeta[position_sec][position_tri][j];

            Z_zetazeta += G[position_tri][i][j] * uv[i] *    Zzetazeta[position_sec][position_tri][j];

         pres_zetazeta += G[position_tri][i][j] * uv[i] *    pressure_3_tri_18[2][position_tri][j];
            
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

    pres_R =      p[0][0] * pres_u + p[0][1] * pres_v + p[0][2] * pres_uu + p[0][3] * pres_uv + p[0][4] * pres_vv + p[0][5] * pres_zetau + p[0][6] * pres_zetav + p[0][7] * pres_zetazeta + p[0][8] * pres_zeta;
    pres_RR =     p[2][0] * pres_u + p[2][1] * pres_v + p[2][2] * pres_uu + p[2][3] * pres_uv + p[2][4] * pres_vv + p[2][5] * pres_zetau + p[2][6] * pres_zetav + p[2][7] * pres_zetazeta + p[2][8] * pres_zeta;
    pres_ZZ =     p[4][0] * pres_u + p[4][1] * pres_v + p[4][2] * pres_uu + p[4][3] * pres_uv + p[4][4] * pres_vv + p[4][5] * pres_zetau + p[4][6] * pres_zetav + p[4][7] * pres_zetazeta + p[4][8] * pres_zeta;
    pres_phiphi = p[7][0] * pres_u + p[7][1] * pres_v + p[7][2] * pres_uu + p[7][3] * pres_uv + p[7][4] * pres_vv + p[7][5] * pres_zetau + p[7][6] * pres_zetav + p[7][7] * pres_zetazeta + p[7][8] * pres_zeta;

    // pres_R =      p[0][0] * pres_u + p[0][1] * pres_v;
    // pres_RR =     p[2][0] * pres_u + p[2][1] * pres_v + p[2][2] * pres_uu + p[2][3] * pres_uv + p[2][4] * pres_vv;
    // pres_ZZ =     p[4][0] * pres_u + p[4][1] * pres_v + p[4][2] * pres_uu + p[4][3] * pres_uv + p[4][4] * pres_vv;
    // pres_phiphi = p[7][0] * pres_u + p[7][1] * pres_v + p[7][2] * pres_uu + p[7][3] * pres_uv + p[7][4] * pres_vv + p[7][5] * pres_zetau + p[7][6] * pres_zetav + p[7][7] * pres_zetazeta + p[7][8] * pres_zeta;

    //R  Z  RR RZ ZZ φR φZ φφ φ
    //0  1  2  3  4  5  6  7  8
    double jacobian =  D * R_;
    double var1 = (pres_R / R_ + pres_RR + pres_ZZ + pres_phiphi / (R_*R_)) * jacobian;

    for (int i = 0; i < dim; i++) {
        delete [] p[i];
    }
    delete [] p;

    return var1;
}






double ***solve_magnetic_field_sparse_LM(double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangles, int num_scalar, double error, time_t time_s) {

    printf("stiffness matrix (3D) assembling...\n");
    fflush(stdout);

    int num_edges = points_sequence[0][9];
    int num_vertices = points_sequence[0][10];
    // printf("num_edges = %d\n", num_edges);

    int num_modes_scalar = num_modes*num_scalar;
    int num_total_DoFs = num_vertices*6*num_modes_scalar;
    int num_boundary_DoFs = num_edges*num_modes_scalar;

    int num_extended_DoFs = num_total_DoFs + num_boundary_DoFs;
 
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

    for (int i = 0; i < num_triangles; i++) {
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
    for (int i = 0; i < num_triangles; i++) {
        // dirichlet boundary condition 
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


    A.setFromTriplets(coeff_A.begin(), coeff_A.end());
    A.makeCompressed();
    vector<T>().swap(coeff_A);

    for (int i = 0; i < num_triangles; i++) {
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
        F_sca_tri_18N[i] = new double *[num_triangles];
        for (int j = 0; j < num_triangles; j++) {
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

double ***solve_pressure_sparse_LM(double ***M_tri_18Nn_18Nn, double **B_tri_18Nn, int **points_sequence, double **boundary_DoFs, int num_modes, int num_triangles, double error, time_t time_s) {

    printf("stiffness matrix (1D) assembling...\n");
    fflush(stdout);

    int num_edges = points_sequence[0][9];
    int num_vertices = points_sequence[0][10];

    int num_total_DoFs = num_vertices*6*num_modes;
    int num_boundary_DoFs = num_edges*num_modes;

    int num_extended_DoFs = num_total_DoFs + num_boundary_DoFs;
    
    SparseMatrix<double> A(num_extended_DoFs, num_extended_DoFs);
    A.setZero();

    VectorXd B(num_extended_DoFs);
    B.setZero();

    VectorXd X(num_extended_DoFs);
    VectorXd X_upper(num_total_DoFs);
    X.setZero();
    X_upper.setZero();

    double **X_re = new double *[1];
    for (int i = 0; i < 1; i++) {
        X_re[i] = new double [num_vertices*6*num_modes]();
    }

    typedef Eigen::Triplet<double> T;
    vector<T> coeff_A;

    Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver1;
    // Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver1;
    // Eigen::LeastSquaresConjugateGradient<SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver1;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver2;
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver3;

    int index_J, index_K, index_j1, index_j2, index_k1;

    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 3; j++) {
            index_J = 6*num_modes*j;
            index_j1 = points_sequence[i][2+3*j]*6*num_modes;
            for (int k = 0; k < 3; k++) {
                index_K = 6*num_modes*k;
                index_k1 = points_sequence[i][2+3*k]*6*num_modes;
                for (int s = 0; s < 6*num_modes; s++) {
                    for (int t = 0; t < 6*num_modes; t++) {
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
    for (int i = 0; i < num_triangles; i++) {
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

                index2 = points_sequence[i][1+3*j]*num_modes;
                index3 = points_sequence[i][2+3*j]*6*num_modes;
                for (int s = 0; s < num_modes; s++) {
                    coeff_A.push_back( T(num_total_DoFs + index2 + s, index3 + s, 1.0) );
                    coeff_A.push_back( T(index3 + s, num_total_DoFs + index2 + s, 1.0) );
                }
                
            }
        }
    }

    A.setFromTriplets(coeff_A.begin(), coeff_A.end());
    A.makeCompressed();
    vector<T>().swap(coeff_A);

    for (int i = 0; i < num_triangles; i++) {
        for (int j = 0; j < 3; j++) {
            index1 = 6*num_modes*j;
            index2 = points_sequence[i][2+3*j]*6*num_modes;
            for (int s = 0; s < 6*num_modes; s++) {
                B(index2 + s) += B_tri_18Nn[i][index1 + s];
            }

            if (points_sequence[i][0+3*j] == 1) {
                index2 = points_sequence[i][1+3*j]*num_modes;
                index3 = points_sequence[i][1+3*j]*num_modes;
                for (int s = 0; s < num_modes; s++) {
                    B(num_total_DoFs + index2 + s) = boundary_DoFs[0][index3 + s];
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
    printf("estimated error = %0.3e\n", abs((A * X).norm() - B.norm()) / (B.norm()));
    timing(time_s);

    // printf(" preconditioning...\n");
    // solver1.compute(A);
    // if (solver1.info() != Eigen::Success) {
    //     printf("solve_magnetic_field_sparse_lagrange_multiplier: decomposition failed...\n");
    //     exit(0);
    // }
    // fflush(stdout); 
    // solver1.setMaxIterations(1999);
    // X = solver1.solve(B);
    // printf("estimated error = %0.3e, iterations = %d\n", solver1.error(), solver1.iterations());
    // timing(time_s);

    for (int i = 0; i < num_total_DoFs; i++) {
        X_upper(i) = X(i);
    }


    printf(" renumbering...\n");
    fflush(stdout);
    for (int i = 0; i < num_total_DoFs; i++) {
        X_re[0][i] = X_upper(i);
    }

    // for (int j = 0; j < 1; j++) {
    //     for (int i = 0; i < num_total_DoFs/1; i++) {
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
    double ***F_sca_tri_18N = new double **[1];
    for (int i = 0; i < 1; i++) {
        F_sca_tri_18N[i] = new double *[num_triangles];
        for (int j = 0; j < num_triangles; j++) {
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




double pressure_initialization(double u, double v, int position_triangle, double **triangles_a_b_c_theta_x0_y0, vector<double> coefficients) {
    
    double theta = triangles_a_b_c_theta_x0_y0[position_triangle][3];
    double cosine = cos(theta);
    double sine = sin(theta);
    double x0 = triangles_a_b_c_theta_x0_y0[position_triangle][4];
    double y0 = triangles_a_b_c_theta_x0_y0[position_triangle][5];
    
    double x = x0 + u * cosine - v * sine;
    double y = y0 + u * sine   + v * cosine;
    double rho = sqrt(x*x + y*y);

    double value = 0;

    for (int i = 0; i < coefficients.size(); i++) {
        value += coefficients[i] * power_nonnegative(rho, i);
    }
    // printf("%0.9f\n", value);
    // exit(0);
    return value;
}
