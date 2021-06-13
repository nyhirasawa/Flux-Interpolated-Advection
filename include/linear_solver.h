#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H
#include "sparse_matrix.h"
#include <vector>

#include "define_float_type.h"

namespace linear_algebra{
//CG法Denseバージョン
void conjugate_gradient(const std::vector<std::vector<MY_FLOAT_TYPE>> &A, const std::vector<MY_FLOAT_TYPE> &b, std::vector<MY_FLOAT_TYPE> &x, int n);
//CG法sparseバージョン
void conjugate_gradient(const sparse_matrix &A, const std::vector<MY_FLOAT_TYPE> &b, std::vector<MY_FLOAT_TYPE> &x, int n, int max_itr, MY_FLOAT_TYPE eps);

//不完全コレスキー分解
void incomplete_cholesky_decomposition(const sparse_matrix &A, sparse_matrix &L, std::vector<MY_FLOAT_TYPE> &d, const int n);
//不完全コレスキー分解によって得られた下三角行列Lと対角行列(の対角成分を保持したベクトルd)から
//(LDL^{T})^{-1}rを計算する関数. 結果はresultに格納される
void calc_LDLt_inv_r(const sparse_matrix& L, const std::vector<MY_FLOAT_TYPE>& d, const std::vector<MY_FLOAT_TYPE>& r, std::vector<MY_FLOAT_TYPE>& result, int n);
//ICCG法sparseバージョン
void incomplete_cholesky_conjugate_gradient(const sparse_matrix &A, const std::vector<MY_FLOAT_TYPE> &b, std::vector<MY_FLOAT_TYPE> &x, int n, int max_itr, MY_FLOAT_TYPE eps);


//gauss_seidel法sparseバージョン
void gauss_seidel(const sparse_matrix_with_diagonal_element &A, const std::vector<MY_FLOAT_TYPE> &b, std::vector<MY_FLOAT_TYPE> &x, int n, int max_itr);
}//namespace linear_algebra
#endif
