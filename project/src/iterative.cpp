#include "slau.h"

void print_big_matr(TYPE **A, int num)
{
    int n = 200 + num;
    TYPE **L;
    TYPE *l;
    l = new TYPE[n];
    L = new TYPE*[n];
    for (int i = 0; i < n; i++)
    {
        L[i] = new TYPE[n];
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
            //   cout << A[1][j] << " ";
                L[i][j] = A[1][j];
            }
            else if (i == j + 1)
            {
            //    cout << A[0][j] << " ";
                L[i][j] = A[0][j];
            }
            else if (i == j - 1)
            {
            //    cout << A[2][j] << " ";
                L[i][j] = A[2][j];
            }
            else
            {
            //    cout << "0" << " ";
                L[i][j] = 0;
            }
        }
        //cout << A[3][i] << endl;
        l[i] = A[3][i];
    }
   // print_sys(L, l, n);
    cout << "I have this :" << endl;
    print_m(L, n);
    cout << endl;
    L = relax_matr_c(L, n, 1);
    cout << "I want to get this : " << endl;
    print_m(L, n);
    delete_m(L, n);
    delete(l);
}

TYPE **create_big_matr(int num)
{
    TYPE **A;
    int n = 200 + num;

    A = new TYPE*[4];
    A[0] = new TYPE[n - 1];
    A[1] = new TYPE[n];
    A[2] = new TYPE[n - 1];
    A[3] = new TYPE[n];
    for (int j = 0; j < n; j++)
    {
        A[1][j] = 4;
        A[0][j] = 1;
        A[2][j] = 1;
        if (j == 0)
            A[3][j] = 6;
        else if (j == n - 1)
            A[3][j] = 9 - 3 * (n % 2);
        else
            A[3][j] = 10 - 2 * ((j + 1) % 2);
    }
    return (A);
}

TYPE **relax_matr_c(TYPE **A, const int size, const TYPE w)
{
    TYPE *x = new TYPE[size];
    TYPE **C = new TYPE*[size];
    TYPE sum;

    for (int i = 0; i < size; i++)
    {
        C[i] = new TYPE[size];
        x[i] = 0;
    }
    for (int h = 0; h < size; h++)
    {
        x[h] = 1;
        if (h != 0)
            x[h - 1] = 0;
        for (int i = 0; i < size; i++)
		{
            C[i][h] = 0;
            sum = 0;
			for (int j = 0; j < size; j++)
            {
                if (j < i)
                    sum += A[i][j] * C[j][h];
                else if(i != j)
				    sum += A[i][j] * x[j];
            }
			C[i][h] = (1 - w) * x[i] - w * sum / A[i][i];
        }
    }
    delete[]x;
    return (C);
}

TYPE **yacoby_matr_c(TYPE **A, const int size)
{
    TYPE **C = new TYPE *[size];
    for (int i = 0; i < size; i++)
    {
        C[i] = new TYPE[size];
        for (int j = 0; j < size; j++)
        {
            if (i != j)
                C[i][j] = -A[i][j] / A[i][i];
            else
                C[i][i] = 0;
        }
    }
    return (C);
}

TYPE **simple_iter_matr_c(TYPE **A, const int size, const TYPE tau)
{
        TYPE **C = new TYPE *[size];
        for (int i = 0; i < size; i++)
        {
            C[i] = new TYPE[size];
            for (int j = 0; j < size; j++)
                C[i][j] = -tau * A[i][j];
            C[i][i] += 1;
        }
    return (C);
}

TYPE *yacoby(TYPE **A, TYPE *b, TYPE *x_real, const int size)
{
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE **C;

    for (int i = 0; i < size; i++)
        x[i] = b[i];
    C = yacoby_matr_c(A, size);
    TYPE norm = norm_1_m(C, size);
    TYPE eps1 = abs(1 - norm) / norm * eps;
    TYPE nev;
    int h = 0;
    while (norm_1_v(diff_v(x1, x, size), size) > eps1)
    {
        for (int i = 0; i < size; i++)
        {
            x1[i] = 0;
            for (int j = 0; j < size; j++)
            {
                if (i != j)
                    x1[i] += -A[i][j] / A[i][i] * x[j];
            }
            x1[i] += b[i] / A[i][i];
        }
        swap_v(&x1, &x);
        if (h == 0)
            nev = norm_1_v(diff_v(x1, x_real, size), size);
        h++;
    }
    TYPE norm_error = norm_1_v(diff_v(x, x_real, size), size);
    cout << setprecision(3) << "|" << setw(13) << scientific << norm_1_v(diff_v(x1, x, size), size)
         << "|" << setw(15) << norm_error;
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    cout << "|" << setw(13) << h << "|" << setw(15)
         << norm << "|" << setw(12) << (int)(log(eps1 / nev) / log(norm)) << "|" << endl;
    delete_m(C, size);
    delete[] x1;
    return (x);
}

TYPE *zeydel(TYPE **A, TYPE *b, TYPE *x_real, const int size)
{
    TYPE **C;
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE *tmp = new TYPE[size];

    for (int i = 0; i < size; i++)
        x[i] = b[i];
    C = relax_matr_c(A, size, 1);
    TYPE norm = norm_1_m(C, size);
    TYPE eps1 = abs(1 - norm) / norm * eps;
    TYPE nev;
    int h = 0;
    while (norm_1_v(diff_v(x1, x, size), size) > eps1)
    {
        for (int i = 0; i < size; i++)
        {
            x1[i] = 0;
            tmp[i] = 0;
            for (int j = 0; j < size; j++)
            {
                if (j < i)
                    x1[i] += -A[i][j] / A[i][i] * tmp[j];
                else if (i != j)
                    x1[i] += -A[i][j] / A[i][i] * x[j];
            }
            x1[i] += b[i] / A[i][i];
            tmp[i] = x1[i];
        }
        swap_v(&x1, &x);
        if (h == 0)
            nev = norm_1_v(diff_v(x_real, x1, size), size);
        h++;
    }
    TYPE norm_error = norm_1_v(diff_v(x, x_real, size), size);
    cout << setprecision(3) << "|" << setw(13) << norm_1_v(diff_v(x1, x, size), size)
         << "|" << setw(15) << norm_error;
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    cout << "|" << setw(13) << h << "|" << setw(15)
         << norm << "|" << setw(12) << (int)(log(eps1 / nev) / log(norm)) << "|" << endl;
    delete_m(C, size);
    delete[] x1;
    delete[] tmp;
    return (x);
}

TYPE *relax(TYPE **A, TYPE *b, TYPE *x_real, const int size, const TYPE w)
{
    TYPE **C;
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE *tmp = new TYPE[size];

    for (int i = 0; i < size; i++)
        x[i] = b[i];
    C = relax_matr_c(A, size, w);
    TYPE norm = norm_1_m(C, size);
    TYPE eps1 = abs(1 - norm) / norm * eps;
    TYPE nev;
    int h = 0;
    while (norm_1_v(diff_v(x1, x, size), size) > eps1)
	{
		for (int i = 0; i < size; i++)
		{
            x1[i] = 0;
			for (int j = 0; j < size; j++)
            {
                if (j < i)
                    x1[i] -= w * A[i][j] / A[i][i] * tmp[j];
                else if(i != j)
				    x1[i] -= w * A[i][j] / A[i][i] * x[j];
            }
			x1[i] += w * b[i] / A[i][i] + (1 - w) * x[i];
            tmp[i] = x1[i];
        }
		swap_v(&x1, &x);
        if (h == 0)
            nev = norm_1_v(diff_v(x1, x_real, size), size);
        h++;
    }
    TYPE norm_error = norm_1_v(diff_v(x, x_real, size), size);
    cout << "|" << setw(13) << w << setprecision(3) <<scientific << "|"  << setw(13) << norm_1_v(diff_v(x1, x, size), size) 
    << "|" << setw(15) << norm_error;
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    cout << "|"<< setw(13) << h << "|" << setw(15) 
    << norm << "|" << setw(12) << (int)(log(eps1 / nev) / log(norm)) << "|" << endl;
    delete_m(C, size);
    delete[] x1;
    delete[] tmp;
    return (x);
}

TYPE *simple_iter(TYPE **A, TYPE *b, TYPE *x_real, const int size, const TYPE tau)
{
    TYPE **C;
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];

    int h = 0;
    for (int i = 0; i < size; i++)
        x[i] = b[i];
    C = simple_iter_matr_c(A, size, tau);
    TYPE norm = norm_1_m(C, size);
    TYPE eps1 = abs(1 - norm) / norm * eps;
    TYPE nev;
    while (norm_1_v(diff_v(x1, x, size), size) > eps1)
    {
        for (int i = 0; i < size; i++)
        {
            x1[i] = 0;
            for (int j = 0; j < size; j++)
            {
                if (i != j)
                    x1[i] += -tau * A[i][j] * x[j];
                else
                    x1[i] += (-tau * A[i][j] + 1) * x[j];
            }
            x1[i] += b[i] * tau;
        }
        swap_v(&x1, &x);
        if (h == 0)
            nev = norm_1_v(diff_v(x1, x_real, size), size);
        h++;
    }
    TYPE norm_error = norm_1_v(diff_v(x, x_real, size), size);
    cout << setprecision(3) << "|" << setw(13) << tau << "|" << setw(13) << norm_1_v(diff_v(x1, x, size), size)
         << "|" << setw(15) << norm_error;
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    cout << "|" << setw(13) << h << "|" << setw(15)
         << norm << "|" << setw(12) << (int)(log(eps1 / nev) / log(norm)) << "|" << endl;
    delete_m(C, size);
    delete[] x1;
    return (x);
}