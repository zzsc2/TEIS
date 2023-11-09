#include "header.hpp"
#include "infrastructure.hpp"

double **read_boundary(const char *file_name);

double **fourier_coefficients_dzeta(double **fourier_coefficients, int nfp);

double **fourier_coefficients_dtheta(double **fourier_coefficients, int nfp);

double ***read_coils_p(const char *file_name, int inflation, double tao);