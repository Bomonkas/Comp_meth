#include "slau.h"

int file_input(TYPE*** A, TYPE** B, string name_file) {
	ifstream file(name_file);
	if (!file) {
		cerr << "Error with opening file\n";
		return (-1);
	}
	string h;
	int size;
    file >> size;
    *A = new TYPE*[size];
    for (int i = 0; i < size; i++)
    {
        (*A)[i] = new TYPE[size];
    }
    *B = new TYPE[size];
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			file >> (*A)[i][j];
		}
		file >> (*B)[i];
	}
	file.close();
    return (size);
}

TYPE norm_1_m(TYPE **A, const int size) {
	if (!A) {
		cout << "Matrix is not exist\n";
		return 0;
	}
	double	norm = 0, sum = 0;;
	for (int i = 0; i < size; i++) {
		sum = 0;
		for (int j = 0; j < size; j++) {
			sum += fabs(A[j][i]);
 		}
		if (sum > norm) {
			norm = sum;
		}
	}
	return norm;
}

TYPE norm_inf_m(TYPE **A, const int size) {
	if (!A) {
		cout << "Matrix is not exist\n";
		return 0;
	}
	double	norm = 0, sum = 0;

	for (int i = 0; i < size; i++) {
		sum = 0;
		for (int j = 0; j < size; j++) {
			sum += fabs(A[i][j]);
		}
		if (sum > norm) {
			norm = sum;
		}
	}
	return norm;
}

void delete_m(TYPE** matr, const int size) {
	for (int i = 0; i < size; i++) {
		delete[] matr[i];
	}
	delete[] matr;
}

void print_sys(TYPE** A, TYPE* B, const int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << fixed << setw(13) << setprecision(4) << A[i][j] << "x[" << j << "]";
		}
		cout << " = " << B[i] << endl;
	}
}

void mult_m_n(TYPE **A, TYPE n, const int size)
{
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			A[i][j] *= n;
}

void print_sys_iter(TYPE** A, TYPE* c, const int size) {
	for (int i = 0; i < size; i++)
	{
		cout << "x[" << i << "] = ";
		for (int j = 0; j < size; j++)
			cout <<  fixed << setw(13) << setprecision(4) << A[i][j] << "x[" << j << "]";
		cout << setw(13) << c[i] << endl;
	}
}

int is_zero(TYPE elem) {
	if (fabs(elem) <= eps)
		return (1);
	return (0);
}

void print_v(TYPE* vec, const int size) {
	if (!vec)
	{
		cout << "Matrix is degenerated" << endl;
		return ;
	}
	for (int i = 0; i < size; i++) {
		cout << vec[i] << " ";
	}
	cout << endl;
}

TYPE *diff_v(TYPE *x, TYPE *y, const int size) {
	TYPE *diff = new TYPE[size];

	for (int i = 0; i < size; i++) {
		diff[i] = x[i] - y[i];
	}
	return diff;
}

TYPE norm_inf_v(TYPE *vec, const int size) {
	TYPE max = fabs(vec[0]);

	for (int i = 1; i < size; i++) {
		if (max < fabs(vec[i]))
			max = fabs(vec[i]);
	}
	return max;
}

void mult_v_n(TYPE *x, TYPE n, int size)
{
	for (int i = 0; i < size; i++)
		x[i] *= n;
}

void swap_v(TYPE **x, TYPE **y, int size)
{
	TYPE *tmp;

	tmp = *x;
	*x = *y;
	*y = tmp;
}

TYPE norm_1_v(TYPE *vec, const int size) {
	TYPE sum = 0.0;
	for (int i = 0; i < size; i++) {
		sum += fabs(vec[i]);
	}
	return sum;
}

TYPE** transp_m(TYPE** matr, int size) {
	TYPE tmp;
	for (int i = 0; i < size; i++) {
		for (int j = i; j < size; j++) {
			tmp = matr[i][j];
			matr[i][j] = matr[j][i];
			matr[j][i] = tmp;
		}
	}
	return matr;
}

TYPE* mult_m_v(TYPE** matr, TYPE* vec, const int size) {
	TYPE* res_vec = NULL;
	res_vec = new TYPE[size];

	for (int i = 0; i < size; i++) {
		TYPE sum = 0;
		for (int j = 0; j < size; j++) {
			sum += matr[i][j] * vec[j];
		}
		if (is_zero(sum)) {
			res_vec[i] = 0;
		}
		else {
			res_vec[i] = sum;
		}
	}
	return res_vec;
}

TYPE** mult_m_m(TYPE** l_matr, TYPE** r_matr, const int size) {
	TYPE** res_matr = NULL;
	res_matr = new TYPE *[size];
	for (int i = 0; i < size; i++) {
		res_matr[i] = new TYPE[size];
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			TYPE sum = 0;
			for (int k = 0; k < size; k++) {
				sum += l_matr[i][k] * r_matr[k][j];
			}
			if (is_zero(sum)) {
				res_matr[i][j] = 0;
			}
			else {
				res_matr[i][j] = sum;
			}
		}
	}
	return res_matr;
}

void print_m(TYPE** matr, const int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << fixed << setprecision(4) << setw(13) << matr[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;
}

int	find_and_swap(TYPE** A, TYPE* B, int k, const int size) {
	if (!A || !B) {
		cout << "Matrix is not exist\n";
		return 1;
	}
	TYPE* tmp_a;
	TYPE tmp_b;
	int  num_max = 0;
	TYPE max = 0.0;

	for (int i = k; i < size; i++) {
		if (fabs(A[i][k]) > max) {
			max = fabs(A[i][k]);
			num_max = i;
		}
	}
	if (is_zero(max)) {
		return 1;
	}
	tmp_a = A[num_max];
	A[num_max] = A[k];
	A[k] = tmp_a;
	tmp_b = B[num_max];
	B[num_max] = B[k];
	B[k] = tmp_b;

	return 0;
}

TYPE* reverse_move(TYPE** A, TYPE* B, const int size) {
	if (!A || !B) {
		cout << "Matrix is not exist\n";
		return nullptr;
	}
	TYPE* result = new TYPE[size];
	for (int i = size - 1; i >= 0; i--) {
		TYPE sum = 0.0;
		if (i != size - 1)
			for (int j = i + 1; j < size; j++) {
				sum += A[i][j] * result[j];
			}
		if (is_zero(A[i][i])) {
			return nullptr;
		}
		else {
			result[i] = (B[i] - sum) / A[i][i];
		}
		if (is_zero(result[i])) {
			result[i] = 0.0;
		}
	}
	return (result);
}

TYPE* gauss(TYPE** A, TYPE* B, const int size) {
	if (!A || !B) {
		cout << "Matrix is not exist\n";
		return nullptr;
	}
	TYPE koff = 0;
	for (int k = 0; k < size; k++) {
		if (find_and_swap(A, B, k, size)) {
			cout << "The matrix is degenerate" << endl;
			return nullptr;
		}
		for (int j = k + 1; j < size; j++) {
			koff = A[j][k] / A[k][k];
			for (int i = k; i < size; i++) {
				if (i == k)
					A[j][i] = 0.0;
				else
					A[j][i] -= koff * A[k][i];
				if (is_zero(A[i][j]))
					A[i][j] = 0.0;
			}
			B[j] -= koff * B[k];
			if (is_zero(B[j])) {
				B[j] = 0.0;
			}
		}
	}
	return reverse_move(A, B, size);
}