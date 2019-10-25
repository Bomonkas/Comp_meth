#include "slau.h"

using namespace std;

int main()
{
    TYPE **A;
    TYPE *b;
    TYPE tau;

    int size = file_input(&A, &b, "test.TXT");
	if (size <= 0)
		return (-1);
    print_sys(A, b, size);
    cout << "Simple iter : " << endl;
    tau = 0.23;
    cout << "|     tau     | discrepancy |     error     | num of iter |     ||C||     | posteriori |" << endl;
    simple_iter(A, b, size, 0.23);
    simple_iter(A, b, size, 0.1);
    simple_iter(A, b, size, 0.33);
    cout << "Yacoby iter : " << endl; 
    cout << "| discrepancy |     error     | num of iter |     ||C||     | posteriori |" << endl;
    yacoby_iter(A, b, size);
    cout << "Zeydel iter : " << endl;
    cout << "| discrepancy |     error     | num of iter |     ||C||     | posteriori |" << endl;
    zey_iter(A, b, size);
    cout << "Relax iter : " << endl;
    cout << "|      w      | discrepancy |     error     | num of iter |     ||C||     | posteriori |" << endl;
    rel_iter(A, b, size, 1.5);
    rel_iter(A, b, size, 0.3);
    cout << "Gauss : " << endl;
    print_v(gauss(A, b, size), size);
    delete(b);
    delete_m(A, size);
    return (0);
}