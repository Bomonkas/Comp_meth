#pragma once
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>

#define TYPE double

#define eps 0.001

using namespace std;

TYPE    **yacoby_matr_c(TYPE **A, const int size);
TYPE    **simple_iter_matr_c(TYPE **A, const int size, const TYPE tau);
TYPE    **relax_matr_c(TYPE **A, const int size, const TYPE w);
TYPE 	*relax(TYPE **A, TYPE *b, TYPE *x_real, const int size, const TYPE w);
TYPE 	*zeydel(TYPE **A, TYPE *b, TYPE *x_real, const int size);
TYPE    *simple_iter(TYPE **A, TYPE *b, TYPE *x_real, const int size, const TYPE tau);
TYPE    *yacoby(TYPE **A, TYPE *B, TYPE *x_real, const int size);
int	    find_and_swap(TYPE** A, TYPE* B, int k, const int size);
TYPE    *reverse_move(TYPE** A, TYPE* B, const int size);
TYPE    *gauss(TYPE** A, TYPE* B, const int size);
void	mult_v_n(TYPE *x, TYPE n, int size);
void	mult_m_n(TYPE **A, TYPE n, const int size);
void    swap_v(TYPE **x, TYPE **y);
TYPE    norm_inf_v(TYPE *vec, const int size);
TYPE    *diff_v(TYPE *x, TYPE *y, const int size);
TYPE    norm_1_v(TYPE *vec, const int size);
TYPE    norm_1_m(TYPE **A, const int size);
TYPE    norm_inf_m(TYPE **A, const int size);
void    print_sys_iter(TYPE** A, TYPE* c, const int size);
int     file_input(TYPE*** A, TYPE** B, string name_file);
TYPE    **inverse_m(TYPE **A, TYPE **I, const int size);
void    delete_m(TYPE** matr, const int size);
void    print_sys(TYPE** A, TYPE* B, const int size);
int     is_zero(TYPE elem);
void    print_v(TYPE* vec, const int size);
TYPE    **transp_m(TYPE** matr, int size);
TYPE    *mult_m_v(TYPE** matr, TYPE* vec, const int size);
TYPE    **mult_m_m(TYPE** l_matr, TYPE** r_matr, const int size);
void    print_m(TYPE** matr, const int size);
TYPE    **create_big_matr(int num);
void    print_big_matr(TYPE **A, int num);
TYPE    **copy_m(TYPE **matr, const int size);
TYPE    *copy_v(TYPE *vec, const int size);
