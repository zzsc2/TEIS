#pragma once
#include <omp.h>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <numeric>
#include <cmath>
#include <complex>
#include <functional>
#include <ctime>
#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Core>
#include <eigen-3.4.0/Eigen/Sparse>
// #include "infrastructure.hpp"

#define pi 3.14159265358979
#define mu0 0.000001256637061435917
#define bs mu0/(4.0*pi)
#define NT 5

using std::vector;
using std::string;
using std::fstream;
using std::accumulate;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::SparseVector;


