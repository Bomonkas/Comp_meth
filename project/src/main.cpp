#include "slau.h"

using namespace std;

int main()
{
    TYPE **A;
    TYPE *b;

    int size = file_input(&A, &b, "tests/T1.TXT");
	if (size <= 0)
		return (-1);
    print_sys(A, b, size);

    cout << "\nSimple iter : " << endl;
    cout << "|     tau     | discrepancy |     error     | num of iter |     ||C||     | posteriori |" << endl;
    TYPE *sim_vec_1 = simple_iter(A, b, size, 0.23);
    TYPE *sim_vec_2 = simple_iter(A, b, size, 0.1);
    TYPE *sim_vec_3 = simple_iter(A, b, size, 0.33);
    TYPE *sim_vec_4 = simple_iter(A, b, size, 0.24);
    cout << "X_1 = ";
    print_v(sim_vec_1, size);
    cout << "X_2 = ";
    print_v(sim_vec_2, size);
    cout << "X_3 = ";
    print_v(sim_vec_3, size);
    cout << "X_4 = ";
    print_v(sim_vec_4, size);

    cout << "\nYacoby iter : " << endl; 
    cout << "| discrepancy |     error     | num of iter |     ||C||     | posteriori |" << endl;
    TYPE *yac_vec = yacoby(A, b, size);
    cout << "X = ";
    print_v(yac_vec, size);

    cout << "\nZeydel iter : " << endl;
    cout << "| discrepancy |     error     | num of iter |     ||C||     | posteriori |\n";
    TYPE *zey_vec = zeydel(A, b, size);
    cout << "X = ";
    print_v(zey_vec, size);

    cout << "\nRelax iter : " << endl;
    cout << "|      w      | discrepancy |     error     | num of iter |     ||C||     | posteriori |" << endl;
    TYPE *rel_vec_1 = relax(A, b, size, 0.8);
    TYPE *rel_vec_2 = relax(A, b, size, 0.4);
    TYPE *rel_vec_3 = relax(A, b, size, 0.1);
    cout << "X_1 = ";
    print_v(rel_vec_1, size);
    cout << "X_2 = ";
    print_v(rel_vec_2, size);
    cout << "X_3 = ";
    print_v(rel_vec_3, size);

    cout << "Gauss : " << endl;
    print_v(gauss(A, b, size), size);
    delete(b);
    delete_m(A, size);
    return (0);
}