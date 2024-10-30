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
// Функция для проверки, является ли матрица вырожденной
bool is_singular(int n, Matrix matrix) {
	// Преобразуем матрицу в формат Eigen
	Eigen::MatrixXd eigen_matrix(n, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j <n; ++j) {
			eigen_matrix(i, j) = matrix.array[i][j];
		}
	}

	// Вычисляем определитель матрицы с помощью Eigen
	double det = eigen_matrix.determinant();

	// Если определитель близок к нулю, матрица вырожденная
	return std::abs(det) < 1e-6;
}

// Функция для генерации случайной невырожденной матрицы
Matrix generate_non_singular_matrix(int n) {
	Matrix temp(n, n);

	// Используем текущее время в качестве seed для генератора случайных чисел
	srand(time(0));

	// Генерируем случайные значения для матрицы
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			temp.array[i][j] = random_double(0.0, 100.0);
		}
	}

	// Проверяем, является ли матрица вырожденной
	while (is_singular(n,temp)) {
		// Если матрица вырожденная, генерируем новую
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
	cout << "Матрица B: " << endl;
	B.ShowMyMatrix();
	Matrix C = generate_non_singular_matrix(n);
	cout << "Матрица C: " << endl;
	C.ShowMyMatrix();
	cout << endl;
	Matrix InvMat = C.inverseMatrix();
	Matrix BC = B * C;
	Matrix A = InvMat * BC;
	cout << endl;
	cout << "Матрица А: " << endl;
	A.ShowMyMatrix();
	cout << endl << endl;
	cout << "Степенной метод: " << endl;
	cout << setprecision(10) << Power_Method(n, A, B) << endl;




	double* sigma = new double[n];
	for (int i = 0; i < n; i++) {;
		cout << "Введите " << i << "-ю координату веткора сигма:" << endl;
		cin >> sigma[i];
	}
	Matrix SIGMA(n, n, sigma);
	cout << endl << endl << "Обратный степенной метод: " << endl;
	Matrix RES = reverse_power_method(n, A, SIGMA);
	

	
	cout << endl; // обратный работает 
	cout << endl;
}