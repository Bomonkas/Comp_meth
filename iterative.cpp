#include "slau.h"

TYPE *yacoby_iter(TYPE **A, TYPE *b, const int size)
{
    TYPE **B = new TYPE *[size];
    TYPE *c = new TYPE[size];
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];

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
        c[i] = b[i] / A[i][i];
        x[i] = c[i];
    }
    TYPE eps1 = (1 - norm_1_m(B, size)) / norm_1_m(B, size) * eps;
    int h = 0;
	while (norm_1_v(diff_v(x1, x, size), size) > eps1)
	{
		for (int i = 0; i < size; i++)
		{
            x1[i] = B[i][0] * x[0];
			for (int j = 1; j < size; j++)
				x1[i] += B[i][j] * x[j];
			x1[i] += c[i];
        }
		swap_v(&x1, &x, size);
        h++;
	} 
    cout << "Num = " << h << endl;
    delete_m(B, size);
    delete[] c;
    delete[] x1;
    return (x);
}

TYPE *zey_iter(TYPE **A, TYPE *b, const int size)
{
    TYPE **B = new TYPE *[size];
    TYPE *c = new TYPE[size];
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
        c[i] = b[i] / A[i][i];
        x[i] = c[i];
    }
    TYPE eps1 = (1 - norm_1_m(B, size)) / norm_1_m(B, size) * eps;
    int h = 0;
    while (norm_1_v(diff_v(x1, x, size), size) > eps1)
	{
		for (int i = 0; i < size; i++)
		{
            x1[i] = 0;
			for (int j = 0; j < size; j++)
            {
                if (j < i)
                    x1[i] += B[i][j] * tmp[j];
                else
				    x1[i] += B[i][j] * x[j];
            }
			x1[i] += c[i];
            tmp[i] = x1[i];
        }
		swap_v(&x1, &x, size);
        h++;
	} 
    cout << "Num = " << h << endl;
    delete_m(B, size);
    delete[] c;
    delete[] x1;
    delete[] tmp;
    return (x);
}

TYPE *rel_iter(TYPE **A, TYPE *b, const int size)
{
    TYPE **B = new TYPE *[size];
    TYPE *c = new TYPE[size];
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
    TYPE *tmp = new TYPE[size];
    TYPE w = 1.5;

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
        c[i] = b[i] / A[i][i];
        x[i] = c[i];
    }
    TYPE eps1 = (1 - norm_1_m(B, size)) / norm_1_m(B, size) * eps;
    int h = 0;
    while (norm_1_v(diff_v(x1, x, size), size) > eps1)
	{
		for (int i = 0; i < size; i++)
		{
            x1[i] = 0;
			for (int j = 0; j < size; j++)
            {
                if (j < i)
                    x1[i] += B[i][j] * tmp[j];
                else
				    x1[i] += B[i][j] * x[j];
            }
			x1[i] += c[i];
            x1[i] = x1[i] + (w - 1) * (x1[i] - x[i]);
            tmp[i] = x1[i];
        }
		swap_v(&x1, &x, size);
        h++;
	} 
    cout << "Num = " << h << endl;
    delete_m(B, size);
    delete[] c;
    delete[] x1;
    delete[] tmp;
    return (x);
}

TYPE *simple_iter(TYPE **A, TYPE *b, const int size)
{
    TYPE **B = new TYPE *[size];
    TYPE *c = new TYPE[size];
    TYPE *x = new TYPE[size];
    TYPE *x1 = new TYPE[size];
	TYPE tau = 0.23;

    for (int i = 0; i < size; i++)
    {
        B[i] = new TYPE[size];
        for (int j = 0; j < size; j++)
        {
            B[i][j] = -tau * A[i][j];
            if (i == j)
                B[i][j] += 1;        
        }
        c[i] = b[i] * tau;
        x[i] = c[i];
    }
    int h = 0;
    TYPE eps1 = abs((1 - norm_1_m(B, size)) / norm_1_m(B, size) * eps);
	while (norm_1_v(diff_v(x1, x, size), size) > eps1)
	{
		for (int i = 0; i < size; i++)
		{
            x1[i] = 0;
			for (int j = 0; j < size; j++)
				x1[i] += B[i][j] * x[j];
			x1[i] += c[i];
        }
        swap_v(&x1, &x, size);
        h++;
	} 
    cout << "Num = " << h << endl;
    delete_m(B, size);
    delete[] c;
    delete[] x1;
    return (x);
}

