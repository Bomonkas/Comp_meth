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
    print_v(simple_iter(A, b, size, 0.23), size);
    simple_iter(A, b, size, 0.1);
    simple_iter(A, b, size, 0.33);
    cout << "\nYacoby iter : " << endl; 
    cout << "| discrepancy |     error     | num of iter |     ||C||     | posteriori |" << endl;
    print_v(yacoby_iter(A, b, size), size);
    cout << "\nZeydel iter : " << endl;
    cout << "| discrepancy |     error     | num of iter |     ||C||     | posteriori |" << endl;
    print_v(zey_iter(A, b, size), size);
    cout << "\nRelax iter : " << endl;
    cout << "|      w      | discrepancy |     error     | num of iter |     ||C||     | posteriori |" << endl;
    print_v(rel_iter(A, b, size, 1.5), size);
    rel_iter(A, b, size, 0.4);
    cout << "Gauss : " << endl;
    print_v(gauss(A, b, size), size);
    delete(b);
    delete_m(A, size);
    return (0);
}