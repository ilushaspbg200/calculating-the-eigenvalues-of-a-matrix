#pragma once
#include <iostream>
#include <vector>
#include "MatrixClass.h"
#include <cmath>
long double Norm(int n, double** array);
double** NormVec(int n, double** array);
double random_double(double min, double max);
void showVec(int n, double** array);
double** Find_Lambda(int n, Matrix first_mat, Matrix second_mat);
double Power_Method(int n, Matrix A, Matrix Diag);