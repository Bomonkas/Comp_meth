#include "slau.h"

using namespace std;

int main()
{
    TYPE **A;
    TYPE *b;
    int size = file_input(&A, &b, "test.TXT");
	if (size <= 0)
		return (-1);
    print_sys(A, b, size);
    cout << "Yacoby iter : " << endl;  
    print_v(yacoby_iter(A, b, size), size);
    cout << "Simple iter : " << endl;
    print_v(simple_iter(A, b, size), size);
    cout << "Zeydel iter : " << endl;
    print_v(zey_iter(A, b, size), size);
    cout << "Relax iter : " << endl;
    print_v(rel_iter(A, b, size), size);
    cout << "Gauss : " << endl;
    print_v(gauss(A, b, size), size);
    delete(b);
    delete_m(A, size);
    return (0);
}