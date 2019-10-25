#include "slau.h"


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
			C[i][h] = (1 - w) * C[i][h] - w * sum / A[i][i];
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

TYPE *yacoby(TYPE **A, TYPE *b, const int size)
{
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE **C;

    for (int i = 0; i < size; i++)
        x[i] = b[i];
    C = yacoby_matr_c(A, size);
    TYPE norm = norm_1_m(C, size);
    TYPE eps1 = (1 - norm) / norm * eps;
    TYPE nev = norm_1_v(diff_v(x1, x, size), size);
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
        h++;
    }
    cout << setprecision(10) << "|" << setw(13) << norm_1_v(diff_v(x1, x, size), size)
         << "|" << setw(15) << eps1 << "|" << setw(13) << h << "|" << setw(15)
         << norm << "|" << setw(12) << (int)(log(eps * (1 - norm) / nev) / log(norm)) << "|" << endl;
    delete_m(C, size);
    delete[] x1;
    return (x);
}

TYPE *zeydel(TYPE **A, TYPE *b, const int size)
{
    TYPE **C;
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE *tmp = new TYPE[size];

    for (int i = 0; i < size; i++)
        x[i] = b[i];
    C = relax_matr_c(A, size, 1);
    TYPE norm = norm_1_m(C, size);
    TYPE eps1 = (1 - norm) / norm * eps;
    TYPE nev = norm_1_v(diff_v(x1, x, size), size);
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
        h++;
    }
    cout << setprecision(10) << "|" << setw(13) << norm_1_v(diff_v(x1, x, size), size)
         << "|" << setw(15) << eps1 << "|" << setw(13) << h << "|" << setw(15)
         << norm << "|" << setw(12) << (int)(log(eps * (1 - norm) / nev) / log(norm)) << "|" << endl;
    delete_m(C, size);
    delete[] x1;
    delete[] tmp;
    return (x);
}

TYPE *relax(TYPE **A, TYPE *b, const int size, const TYPE w)
{
    TYPE **C;
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE *tmp = new TYPE[size];

    for (int i = 0; i < size; i++)
        x[i] = b[i];
    C = relax_matr_c(A, size, w);
    TYPE norm = norm_1_m(C, size);
    TYPE eps1 = (1 - norm) / norm * eps;
    TYPE nev = norm_1_v(diff_v(x1, x, size), size);
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
        h++;
	} 
    cout << "|" << setw(13) << w  << setprecision(10) << "|"  << setw(13) << norm_1_v(diff_v(x1, x, size), size) 
    << "|" << setw(15) << eps1 << "|"<< setw(13) << h << "|" << setw(15) 
    << norm << "|" << setw(12) << (int)(log(eps * (1 - norm) / nev) / log(norm)) << "|" << endl;
    delete_m(C, size);
    delete[] x1;
    delete[] tmp;
    return (x);
}

TYPE *simple_iter(TYPE **A, TYPE *b, const int size, const TYPE tau)
{
    TYPE **C;
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];

    int h = 0;
    for (int i = 0; i < size; i++)
        x[i] = b[i];
    C = simple_iter_matr_c(A, size, tau);
    TYPE norm = norm_1_m(C, size);
    TYPE eps1 = (1 - norm) / norm * eps;
    TYPE nev = norm_1_v(diff_v(x1, x, size), size);
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
        h++;
    }
    cout << setprecision(10) << "|" << setw(13) << tau << "|" << setw(13) << norm_1_v(diff_v(x1, x, size), size)
         << "|" << setw(15) << eps1 << "|" << setw(13) << h << "|" << setw(15)
         << norm << "|" << setw(12) << (int)(log(eps * (1 - norm) / nev) / log(norm)) << "|" << endl;
    delete_m(C, size);
    delete[] x1;
    return (x);
}