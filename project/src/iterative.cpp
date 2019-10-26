#include "slau.h"

TYPE *yacoby_iter(TYPE **A, TYPE *b, const int size)
{
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE **С = new TYPE *[size];

    for (int i = 0; i < size; i++)
    {
        С[i] = new TYPE[size];
        for (int j = 0; j < size; j++)
        {
            if (i != j)
                С[i][j] = -A[i][j] / A[i][i];
            else
                С[i][i] = 0;
        }
        x[i] = b[i];
    }
    TYPE norm = norm_1_m(С, size);
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
    cout << setprecision(10) << "|"  << setw(13) << norm_1_v(diff_v(x1, x, size), size) 
    << "|" << setw(15) << eps1 << "|"<< setw(13) << h << "|" << setw(15) 
    << norm << "|"<< setw(12) << (int)(log(eps * (1 - norm) / nev) / log(norm)) << "|" << endl;
    delete_m(С, size);
    delete[] x1;
    return (x);
}

TYPE **relax_m(TYPE **A, const int size, const TYPE w)
{
    TYPE **C;
    TYPE **D;
    TYPE **L;

    D = new TYPE*[size];
    L = new TYPE*[size];
    C = new TYPE*[size];
    for (int i = 0; i < size; i++)
    {
        C[i] = new TYPE[size];
        D[i] = new TYPE[size];
        L[i] = new TYPE[size];
    }
}

TYPE *zey_iter(TYPE **A, TYPE *b, const int size)
{
    TYPE **C = new TYPE *[size];
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE *tmp = new TYPE[size];

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
        x[i] = b[i];
    }
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
                else if (i!= j)
                    x1[i] += -A[i][j] / A[i][i] * x[j];
            }
			x1[i] += b[i] / A[i][i];
            tmp[i] = x1[i];
        }
		swap_v(&x1, &x);
        h++;
	}
    cout << setprecision(10) << "|"  << setw(13) << norm_1_v(diff_v(x1, x, size), size) 
    << "|" << setw(15) << eps1 << "|"<< setw(13) << h << "|" << setw(15) 
    << norm << "|"<< setw(12) << (int)(log(eps * (1 - norm) / nev) / log(norm)) << "|" << endl;
    delete_m(C, size);
    delete[] x1;
    delete[] tmp;
    return (x);
}

TYPE *rel_iter(TYPE **A, TYPE *b, const int size, const TYPE w)
{
    TYPE **C = new TYPE *[size];
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE *tmp = new TYPE[size];

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
        x[i] = b[i];
    }
    TYPE norm = norm_1_m(C, size);
    TYPE eps1 = (1 - norm) / norm * eps;
    TYPE nev = norm_1_v(diff_v(x1, x, size), size);
    int h = 0;
    while (norm_1_v(diff_v(x1, x, size), size) > eps1)
	{
		for (int i = 0; i < size; i++)
		{
            x1[i] = 0;
			for (int j = 0; j < i; j++)
            {
				x1[i] -= w * A[i][j] / A[i][i] * x1[j];
            }
            x1[i] += (1 - w) * x[i] + w * b[i] / A[i][i];
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
    TYPE **C = new TYPE *[size];
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];

    for (int i = 0; i < size; i++)
    {
        C[i] = new TYPE[size];
        for (int j = 0; j < size; j++)
        {
            C[i][j] = -tau * A[i][j];     
        }
        C[i][i] += 1;
        x[i] = b[i];
    }

    int h = 0;
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
                    x1[i] += -tau * A[i][j] * x[j] + 1;
            }
			x1[i] += b[i] * tau;
        }
        swap_v(&x1, &x);
        h++;
	} 
    cout << setprecision(10) << "|" << setw(13) << tau  << "|"  << setw(13) << norm_1_v(diff_v(x1, x, size), size) 
    << "|" << setw(15) << eps1 << "|"<< setw(13) << h << "|" << setw(15) 
    << norm << "|"<< setw(12) << (int)(log(eps * (1 - norm) / nev) / log(norm)) << "|" << endl;
    delete_m(C, size);
    delete[] x1;
    return (x);
}

