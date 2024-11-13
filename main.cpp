#include <iostream>
#include <vector>
#include "MatrixClass.h"
#include "PowerMethod.h"
#include "reverse_power_method.h"
#include "QR_method.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <Eigen/Dense>

//double random_double(double min, double max) {
//	return min + static_cast<double>(rand()) / (RAND_MAX / (max - min));
//}
Matrix DiagMatr(int n) {
	Matrix temp(n, n);
	for (int i = 0; i < n; i++) {
		temp.array[i][i] = random_double(0.0, 100.0);
	}
	return temp;
}

bool is_singular(int n, Matrix matrix) {

	Eigen::MatrixXd eigen_matrix(n, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j <n; ++j) {
			eigen_matrix(i, j) = matrix.array[i][j];
		}
	}

	
	double det = eigen_matrix.determinant();

	
	return std::abs(det) < 1e-6;
}


Matrix generate_non_singular_matrix(int n) {
	Matrix temp(n, n);
	srand(time(0));

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			temp.array[i][j] = random_double(0.0, 100.0);
		}
	}


	while (is_singular(n,temp)) {

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				temp.array[i][j] = random_double(0.0, 100.0);
			}
		}
	}

	return temp;
}
using namespace std;
int main() {
	setlocale(LC_ALL, "rus");
	srand(time(0));
	int n = 5;
	Matrix B = DiagMatr(n);
	cout << "Matrix B: " << endl;
	B.ShowMyMatrix();
	Matrix C = generate_non_singular_matrix(n);
	cout << "Matrix C: " << endl;
	C.ShowMyMatrix();
	cout << endl;
	Matrix InvMat = C.inverseMatrix();
	Matrix BC = B * C;
	Matrix A = InvMat * BC;
	cout << endl;
	cout << "Matrix À: " << endl;
	A.ShowMyMatrix();
	cout << endl << endl;
	cout << "Power method: " << endl;
	cout << setprecision(10) << Power_Method(n, A, B) << endl;




	double* sigma = new double[n];
	for (int i = 0; i < n; i++) {;
		cout << "Enter " << i << "the ith coordinate of the sigma vector:" << endl;
		cin >> sigma[i];
	}
	Matrix SIGMA(n, n, sigma);
	cout << endl << endl << "Reverse power method: " << endl;
	Matrix RES = reverse_power_method(n, A, SIGMA);
	

	
	cout << endl;  

}
