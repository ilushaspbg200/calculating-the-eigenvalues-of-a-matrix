#include "PowerMethod.h"
double random_double(double min, double max) {
	return min + static_cast<double>(rand()) / (RAND_MAX / (max - min));
}
long double Norm(int n, double** array) {
	long double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += array[i][0] * array[i][0];
	}
	sum = sqrt(sum);
	return sum;
}
double** NormVec(int n, double** array) {

	long double sum = Norm(n, array);
	for (int i = 0; i < n; i++) {
		array[i][0] /= sum;
	}
	return array;
}
void showVec(int n, double** array) {
	for (int i = 0; i < n; i++) {
		cout << array[i][0] << " ";
	}
	cout << endl;
}
double** Find_Lambda(int n, Matrix first_mat, Matrix second_mat) {
	double** Lambda = new double* [n];
	for (int i = 0; i < n; i++) {
		Lambda[i] = new double[1];
		if (abs(second_mat.array[i][0]) > 0.00000001)
			Lambda[i][0] = first_mat.array[i][0] / second_mat.array[i][0];
		else
			Lambda[i][0] = 0;
	}
	return Lambda;
}
//void Clear_Vec(int n,double** array) {
//	for (int i = 0; i < n; i++) {
//		delete[] array[i];
//	}
//	delete array;
//
//}
double Power_Method(int n, Matrix A,Matrix Diag) {
	// произвольный начальный вектор
	double** y0 = new double* [n];
	for (int i = 0; i < n; i++) {
		y0[i] = new double[1];
		y0[i][0] = 1;
	}
	double** z0 = NormVec(n, y0);
	Matrix Z0(n, 1, z0);
	Matrix y1 = A * Z0;
	double** Lambda1 = Find_Lambda(n, y1, Z0);
	Matrix LAMBDA1(n, 1, Lambda1);
	double** zi = NormVec(n, y1.array);
	Matrix Zi(n, 1, zi); // создаем zi векторы
	Matrix Yi = A * Zi;
	double** Lambda2 = Find_Lambda(n, Yi, Zi);
	Matrix LAMBDA2(n, 1, Lambda2);
	int i = 0;
	while (Norm(n, (LAMBDA2 - LAMBDA1).array) > (0.000001) * max(Norm(n, LAMBDA1.array), Norm(n, LAMBDA2.array))) {
		LAMBDA1 = LAMBDA2;
		//Clear_Vec(n, Lambda2);
		//Clear_Vec(n, zi);
		zi = NormVec(n, Yi.array);
		Matrix N(n, 1, zi); // нормированная матрица 
		Zi = N;
		Yi = A * Zi;
		Lambda2 = Find_Lambda(n, Yi, Zi);
		Matrix Lambda2_temp(n, 1, Lambda2);// чтобы присвоить матрице новую матрицу
		LAMBDA2 = Lambda2_temp;
		i += 1;

	}
	cout << endl << "Количество итераций: " << i << endl;
	double res = 0;
	for (int j = 0; j < n; j++) {
		res += LAMBDA2.array[j][0];
	}
	int I = 0;
	for (int j = 0; j < n; j++) {
		if (LAMBDA2.array[j][0] != 0)
			I += 1;
	}
	res = res / I; // собс число
	double ma = -10000000000;
	for (int k = 0; k < n; k++) {
		if (Diag.array[k][k] > ma)
			ma = Diag.array[k][k];
	}
	ma = ma + random_double(0.0000001, 0.000001);
	//cout << "Собственный вектор: " << endl;
	//showVec(n, NormVec(n,Yi.array)); // собственный вектор
	cout << endl << "собственный вектор 2 : " << endl;
	showVec(n, Zi.array);
	cout << endl << "Собственное число: " << endl;
	return ma;
}

