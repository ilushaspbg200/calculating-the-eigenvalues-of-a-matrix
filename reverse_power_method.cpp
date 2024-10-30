#include "reverse_power_method.h"

Matrix generate_one_mat(int n) {
	Matrix A(n, n);
	for (int i = 0; i < n; i++) {
		A.array[i][i] = 1;
	}
	return A;
}
Matrix generate_Mew(int n, Matrix Z, Matrix Y) {
	Matrix temp(n, 1);
	for (int i = 0; i < n; i++) {
		if (abs(Y.array[i][0]) > 0.00000001)
			temp.array[i][0] = Z.array[i][0] / Y.array[i][0];
		else
			temp.array[i][0] = 0;
	}
	return temp;
}

int Find_I(int n, Matrix Mew) {
	int I = 0;
	for (int i = 0; i < n; i++) {
		if (Mew.array[i][0] != 0)
			I += 1;
	}
	return I;
}
double Sum_Mew(int n, Matrix Mew) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += Mew.array[i][0];
	}
	sum = sum / Find_I(n, Mew);
	return sum;
}

Matrix gauss(int n,Matrix A, Matrix b) {
	for (int i = 0; i < n; ++i) {
		// ѕоиск максимального элемента
		int max_row = i;
		for (int k = i + 1; k < n; ++k) {
			if (std::abs(A.array[k][i]) > std::abs(A.array[max_row][i])) {
				max_row = k;
			}
		}
		std::swap(A.array[i], A.array[max_row]);
		std::swap(b.array[i], b.array[max_row]);
		for (int k = i + 1; k < n; ++k) {
			double factor = A.array[k][i] / A.array[i][i];
			for (int j = i; j < n; ++j) {
				A.array[k][j] -= factor * A.array[i][j];
			}
			b.array[k][0] -= factor * b.array[i][0];
		}
	}
	Matrix x(n, 1);
	for (int i = n - 1; i >= 0; --i) {
		x.array[i][0] = b.array[i][0];
		for (int j = i + 1; j < n; ++j) {
			x.array[i][0] -= A.array[i][j] * x.array[j][0];
		}
		x.array[i][0] /= A.array[i][i];
	}
	return x;
}

Matrix reverse_power_method(int n,Matrix A, Matrix sigma0)
{
	Matrix Vec(n, n); // здесь собственные векторы
	for (int i = 0; i < n; i++) {
		double sigmai_1 = sigma0.array[i][0];
		// произвольный начальный вектор
		double** y0 = new double* [n];
		for (int j = 0; j < n; j++) {
			y0[j] = new double[1];
			y0[j][0] = 1;
		}
		double** z0 = NormVec(n, y0);
		Matrix Z0(n, 1, z0);
		Matrix Yi = gauss(n, (A - generate_one_mat(n)* sigmai_1), Z0); // проверить работает ли умножение на число 
		double** z1 = NormVec(n, Yi.array);
		Matrix Z1(n, 1, z1);
		Matrix Mewi = generate_Mew(n, Z0, Yi);

		double sigmai = sigmai_1 + Sum_Mew(n, Mewi);
		while (abs(sigmai - sigmai_1) > 0.000001) {
			sigmai_1 = sigmai;
			Yi = gauss(n, (A - generate_one_mat(n) * sigmai), Z1);
			Mewi = generate_Mew(n, Z1, Yi);
			z1 = NormVec(n, Yi.array);
			Matrix N(n, 1, z1);
			Z1 = N;
			sigmai = sigmai_1 + Sum_Mew(n, Mewi);
		}
		sigma0.array[i][0] = sigmai;
		for (int j = 0; j < n; j++) {
			Vec.array[j][i] = Z1.array[j][0];
		}
	}
	cout << endl << "—обственные векторы: " << endl;
	Vec.ShowMyMatrix();
	cout << endl;
	cout << "—обственные значени€: " << endl;
	for (int i = 0; i < n; i++) {
		cout << sigma0.array[i][0] << " ";
	}
	cout << endl;
	return sigma0;
}