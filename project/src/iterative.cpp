#include "slau.h"

TYPE *yacoby_iter(TYPE **A, TYPE *b, const int size)
{
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE **B = new TYPE *[size];

    for (int i = 0; i < size; i++)
    {
        B[i] = new TYPE[size];
        for (int j = 0; j < size; j++)
        {
            if (i != j)
                B[i][j] = -A[i][j] / A[i][i];
            else
                B[i][i] = 0;
        }
        x[i] = b[i];
    }
    TYPE norm = norm_1_m(B, size);
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
                if (i !=j )
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
    delete_m(B, size);
    delete[] x1;
    return (x);
}

TYPE *zey_iter(TYPE **A, TYPE *b, const int size)
{
    TYPE **B = new TYPE *[size];
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE *tmp = new TYPE[size];

    for (int i = 0; i < size; i++)
    {
        B[i] = new TYPE[size];
        for (int j = 0; j < size; j++)
        {
            if (i != j)
                B[i][j] = -A[i][j] / A[i][i];
            else
                B[i][i] = 0;
        }
        x[i] = b[i];
    }
    TYPE norm = norm_1_m(B, size);
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
    delete_m(B, size);
    delete[] x1;
    delete[] tmp;
    return (x);
}

TYPE *rel_iter(TYPE **A, TYPE *b, const int size, const TYPE w)
{
    TYPE **B = new TYPE *[size];
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE *tmp = new TYPE[size];

    for (int i = 0; i < size; i++)
    {
        B[i] = new TYPE[size];
        for (int j = 0; j < size; j++)
        {
            if (i != j)
                B[i][j] = -A[i][j] / A[i][i];
            else
                B[i][i] = 0;
        }
        x[i] = b[i];
    }
    TYPE norm = norm_1_m(B, size);
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
                    x1[i] += -A[i][j] / A[i][i] * tmp[j];
                else if(i != j)
				    x1[i] += -A[i][j] / A[i][i] * x[j];
            }
			x1[i] += b[i] / A[i][i];
            x1[i] = x1[i] + (w - 1) * (x1[i] - x[i]);
            tmp[i] = x1[i];
        }
		swap_v(&x1, &x);
        h++;
	} 
    cout << "|" << setw(13) << w  << setprecision(10) << "|"  << setw(13) << norm_1_v(diff_v(x1, x, size), size) 
    << "|" << setw(15) << eps1 << "|"<< setw(13) << h << "|" << setw(15) 
    << norm << "|" << setw(12) << (int)(log(eps * (1 - norm) / nev) / log(norm)) << "|" << endl;
    delete_m(B, size);
    delete[] x1;
    delete[] tmp;
    return (x);
}

TYPE *simple_iter(TYPE **A, TYPE *b, const int size, const TYPE tau)
{
    TYPE **B = new TYPE *[size];
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];

    for (int i = 0; i < size; i++)
    {
        B[i] = new TYPE[size];
        for (int j = 0; j < size; j++)
        {
            B[i][j] = -tau * A[i][j];
            if (i == j)
                B[i][j] += 1;        
        }
        x[i] = b[i];
    }
    int h = 0;
    TYPE norm = norm_1_m(B, size);
    TYPE eps1 = (1 - norm) / norm * eps;
    TYPE nev = norm_1_v(diff_v(x1, x, size), size);
	while (norm_1_v(diff_v(x1, x, size), size) > eps1)
	{
		for (int i = 0; i < size; i++)
		{
            x1[i] = 0;
			for (int j = 0; j < size; j++)
				x1[i] += B[i][j] * x[j];
			x1[i] += b[i] * tau;
        }
        swap_v(&x1, &x);
        h++;
	} 
    cout << setprecision(10) << "|" << setw(13) << tau  << "|"  << setw(13) << norm_1_v(diff_v(x1, x, size), size) 
    << "|" << setw(15) << eps1 << "|"<< setw(13) << h << "|" << setw(15) 
    << norm_1_m(B, size) << "|"<< setw(12) << (int)(log(eps * (1 - norm) / nev) / log(norm)) << "|" << endl;
    delete_m(B, size);
    delete[] x1;
    return (x);
}

