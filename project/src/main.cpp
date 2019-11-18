#include "slau.h"

using namespace std;

int main()
{
    TYPE **A;
    TYPE *b;

    int size = file_input(&A, &b, "tests/ITER1.TXT");
	if (size <= 0)
		return (-1);
    
    TYPE *X = gauss(A, b, size);
    cout << "gauss : ";
    print_v(X, size);

    TYPE *x_real = new TYPE[size];
    x_real = X;
    print_sys(A, b, size);
    // for (int i = 0; i < size; i++) {
    //     cout << "x_real[" << i << "] = ";
    //     cin >> x_real[i];
    // }
    cout << "x_real = ";
    print_v(x_real, size);

    // cout << "\nSimple iter : " << endl;
    // cout << "|     tau     | discrepancy |     error     | num of iter |     ||C||     |   apriory  |" << endl;
    // TYPE *sim_vec_1 = simple_iter(A, b, x_real, size, 1.0);
    // TYPE *sim_vec_2 = simple_iter(A, b, x_real, size, 0.1);
    // TYPE *sim_vec_3 = simple_iter(A, b, x_real, size, 0.01);
    // TYPE *sim_vec_4 = simple_iter(A, b, x_real, size, 0.001);
    // cout << "X_1 = ";
    // print_v(sim_vec_1, size);
    // cout << "X_2 = ";
    // print_v(sim_vec_2, size);
    // cout << "X_3 = ";
    // print_v(sim_vec_3, size);
    // cout << "X_4 = ";
    // print_v(sim_vec_4, size);

    // cout << "\nYacoby iter : " << endl; 
    // cout << "| discrepancy |     error     | num of iter |     ||C||     |   apriory  |" << endl;
    // TYPE *yac_vec = yacoby(A, b, x_real, size);
    // cout << "X = ";
    // print_v(yac_vec, size);

    // cout << "\nZeydel iter : " << endl;
    // cout << "| discrepancy |     error     | num of iter |     ||C||     |   apriory  |\n";
    // TYPE *zey_vec = zeydel(A, b, x_real, size);
    // cout << "X = ";
    // print_v(zey_vec, size);

    cout << "\nRelax iter : " << endl;
    cout << "|      w      | discrepancy |     error     | num of iter |     ||C||     |   apriory  |" << endl;
    TYPE *rel_vec_1 = relax(A, b, x_real, size, 1);
    TYPE *rel_vec_2 = relax(A, b, x_real,  size, 0.4);
    TYPE *rel_vec_3 = relax(A, b, x_real, size, 0.8);
    TYPE *rel_vec_4 = relax(A, b, x_real, size, 1.2);
    TYPE *rel_vec_5 = relax(A, b, x_real, size, 1.6);
    cout << "X_1 = ";
    print_v(rel_vec_1, size);
    cout << "X_2 = ";
    print_v(rel_vec_2, size);
    cout << "X_3 = ";
    print_v(rel_vec_3, size);
    cout << "X_4 = ";
    print_v(rel_vec_4, size);
    cout << "X_5 = ";
    print_v(rel_vec_5, size);

    print_big_matr(create_big_matr(-195), -195);
    delete(b);
    delete_m(A, size);
    return (0);
}