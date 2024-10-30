#pragma once
#include <iostream>
#include <vector>
#include "MatrixClass.h"
#include <cmath>
#include <iomanip>
#include "PowerMethod.h"

Matrix reverse_power_method(int n, Matrix A, Matrix sigma0);
Matrix gauss(int n, Matrix A, Matrix b);
Matrix generate_one_mat(int n);
