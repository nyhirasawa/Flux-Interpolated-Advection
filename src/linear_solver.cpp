#include "linear_solver.h"
#include <iostream>
#include <fstream>
#include "define_float_type.h"

namespace linear_algebra{
//CG法Denseバージョン
void conjugate_gradient(const std::vector<std::vector<MY_FLOAT_TYPE>> &A, const std::vector<MY_FLOAT_TYPE> &b, std::vector<MY_FLOAT_TYPE> &x, int n){
    std::vector<MY_FLOAT_TYPE> error(n), direction(n), Ap(n);
    for(int i=0;i<n; i++){
        x[i]=0.0;
    }
    //初期error
    for(int ix=0;ix<n;ix++){
        MY_FLOAT_TYPE Ax_tmp=0.0;
        for(int iy=0;iy<n;iy++){
            Ax_tmp+=A[ix][iy]*x[iy];
        }
        error[ix]=b[ix]-Ax_tmp;
        direction[ix]=error[ix];
    }
    //メインの計算部分
    for(int i=0;i<100000;i++){
        for(int ix=0;ix<n;ix++){
            Ap[ix]=0.0;
            for(int iy=0;iy<n;iy++){
                Ap[ix]+=A[ix][iy]*direction[iy];
            }
        }
        //修正係数coeffの計算
        MY_FLOAT_TYPE errerr_prev=0.0, diry=0.0;
        for(int ix=0;ix<n;ix++){
            errerr_prev+=error[ix]*error[ix];
            diry+=Ap[ix]*direction[ix];
        }
        MY_FLOAT_TYPE coeff_0=errerr_prev/diry;
        for(int ix=0;ix<n;ix++){
            //解の近似値を計算
            x[ix]+=coeff_0*direction[ix];
            //errorの計算
            error[ix]-=coeff_0*Ap[ix];
        }
        //directionの計算
        MY_FLOAT_TYPE errerr_next=0.0;
        for(int ix=0;ix<n;ix++){
            errerr_next+=error[ix]*error[ix];
        }
        //終了条件
        if(errerr_next<0.01*0.01){
            return;
        }
        //directionの修正係数coeff_0の計算
        MY_FLOAT_TYPE coeff_dir=errerr_next/errerr_prev;
        for(int ix=0;ix<n;ix++){
            direction[ix]=error[ix]+coeff_dir*direction[ix];
        }
    }
}

//CG法sparseバージョン
void conjugate_gradient(const sparse_matrix &A, const std::vector<MY_FLOAT_TYPE> &b, std::vector<MY_FLOAT_TYPE> &x, int n, int max_itr, MY_FLOAT_TYPE eps){
    std::vector<MY_FLOAT_TYPE> error(n), direction(n), Adir(n);

    for(int i=0;i<n; i++){
        x[i]=0.0;
    }
    //初期error
    std::vector<MY_FLOAT_TYPE> Ax_tmp(n);
    mat_vec_product(A, x, Ax_tmp, n);
    for(int ix=0;ix<n;ix++){
        error[ix]=b[ix]-Ax_tmp[ix];
        direction[ix]=error[ix];
    }
    //初期エラーの大きさ
    MY_FLOAT_TYPE errerr_0=0.0;
    for(int ix=0;ix<n;ix++){
        errerr_0+=error[ix]*error[ix];
    }
    //メインの計算部分
    for(int i=0;i<max_itr;i++){
        mat_vec_product(A, direction, Adir, n);

        //解の修正係数coeff_0の計算
        MY_FLOAT_TYPE errerr_prev=0.0, diry=0.0;
        for(int ix=0;ix<n;ix++){
            errerr_prev+=error[ix]*error[ix];
            diry+=Adir[ix]*direction[ix];
        }
        MY_FLOAT_TYPE coeff_0=errerr_prev/diry;
        for(int ix=0;ix<n;ix++){
            //解の近似値を計算
            x[ix]+=coeff_0*direction[ix];
            //errorの計算
            error[ix]-=coeff_0*Adir[ix];
        }
        //directionの計算
        MY_FLOAT_TYPE errerr_next=0.0;
        for(int ix=0;ix<n;ix++){
            errerr_next+=error[ix]*error[ix];
        }
        //終了条件
        if((errerr_next/errerr_0)<eps*eps){
/*
            std::ofstream writing_file;
            writing_file.open("iteration_num_CG_512_2d.dat", std::ios::app);
            writing_file << i << std::endl;
            writing_file.close();
*/
            return;
        }
/*
        if(errerr_next<eps*eps){
            return;
        }
*/
        //directionの修正係数coeff_0の計算
        MY_FLOAT_TYPE coeff_dir=errerr_next/errerr_prev;
        for(int ix=0;ix<n;ix++){
            direction[ix]=error[ix]+coeff_dir*direction[ix];
        }
    }
/*
    std::ofstream writing_file;
    writing_file.open("iteration_num_CG_512_2d.dat", std::ios::app);
    writing_file << max_itr << std::endl;
    writing_file.close();
*/
}

//不完全コレスキー分解
void incomplete_cholesky_decomposition(const sparse_matrix &A, sparse_matrix &L, std::vector<MY_FLOAT_TYPE> &d, const int n){
    int irow=0;
    int ielem_A=0;
    int ielem_L=0;
    d.resize(n);
    //全ての行を走査するループ
    for(auto itr=A.cumulative_num_nonzero_element_in_row.begin()+1;itr!=A.cumulative_num_nonzero_element_in_row.end();itr++){
        //各行における非ゼロな値を持つ列を舐めるループ
        for(int ic=0;ic<(*(itr)-*(itr-1));ic++){
            //行列Lは下三角行列なので上三角成分は0
            if(A.column_index[ielem_A] > irow){
                ielem_A++;
                continue;
            }
            else{
                //L[irow, A.column_index[ielem_A]]の計算
                L.input_element(irow,A.column_index[ielem_A],A.element_value[ielem_A]);
                //irow行目を走るループ
                for(int k0=L.cumulative_num_nonzero_element_in_row[irow];k0<L.cumulative_num_nonzero_element_in_row[irow+1]-1;k0++){
                    //A.column_index[ielem_A]行目を走るループ
                    for(int k1=L.cumulative_num_nonzero_element_in_row[(A.column_index[ielem_A])];k1<L.cumulative_num_nonzero_element_in_row[(A.column_index[ielem_A])+1]-1;k1++){
                        if(L.column_index[k0]==L.column_index[k1]){
                            L.element_value[ielem_L]-=L.element_value[k0]*d[L.column_index[k0]]*L.element_value[k1];
                        }
                    }
                }
                //対角成分の逆数がd[i]
                if(irow==A.column_index[ielem_A]){
                    d[irow]=1.0/L.element_value[ielem_L];
                }
                ielem_A++;
                ielem_L++;
            }
        }
        irow++;
    }
}

//不完全コレスキー分解によって得られた下三角行列Lと対角行列(の対角成分を保持したベクトルd)から
//(LDL^{T})^{-1}rを計算する関数. 結果はresultに格納される
void calc_LDLt_inv_r(const sparse_matrix& L, const std::vector<MY_FLOAT_TYPE>& d, const std::vector<MY_FLOAT_TYPE>& r, std::vector<MY_FLOAT_TYPE>& result, int n){
    result.resize(n);
    std::vector<MY_FLOAT_TYPE> DLTres(n);
    //DL^{T} result を求める(つまり, L*DLTres=r をDLTresについてとく)
    //全ての行を走査するループ
    int irow=0;
    for(auto itr=L.cumulative_num_nonzero_element_in_row.begin()+1;itr!=L.cumulative_num_nonzero_element_in_row.end();itr++){
        //各行における非ゼロな値を持つ列を舐めるループ
        MY_FLOAT_TYPE LDLTres=r[irow];
        for(int ielem_L=*(itr-1);ielem_L<*(itr)-1;ielem_L++){
            LDLTres-=L.element_value[ielem_L]*DLTres[L.column_index[ielem_L]];
        }
        //対角成分で割る
        DLTres[irow]=LDLTres/L.element_value[*(itr)-1];
        irow++;
    }

    //result を求める(つまり, DL^{T}*result=DLTres をresultについてとく)
    for(int i=0;i<n;i++){
        result[i]=DLTres[i];
    }
    //全ての列を逆向きに走査するループ
    int icol=n-1;
    for(auto itr=L.cumulative_num_nonzero_element_in_row.end()-1;itr!=L.cumulative_num_nonzero_element_in_row.begin();itr--){
        //各列における非ゼロな値を持つ行を舐めるループ
        for(int ielem_L=*(itr)-1;ielem_L>=*(itr-1);ielem_L--){
            //Lは転置しているので行のindexが列に, 列のindexが行のindexになる
            irow=L.column_index[ielem_L];
            if(irow!=icol){
                result[irow]-=L.element_value[ielem_L]*d[irow]*result[icol];
            }
        }
        //対角成分で割る
        icol--;
    }
}

//ICCG法sparseバージョン
void incomplete_cholesky_conjugate_gradient(const sparse_matrix &A, const std::vector<MY_FLOAT_TYPE> &b, std::vector<MY_FLOAT_TYPE> &x, int n, int max_itr, MY_FLOAT_TYPE eps){
    //前処理として不完全コレスキー分解を行う
    sparse_matrix L(n, n);
    std::vector<MY_FLOAT_TYPE> d(n);
    incomplete_cholesky_decomposition(A, L, d, n);


    std::vector<MY_FLOAT_TYPE> error(n), direction(n), Adir(n);
    for(int i=0;i<n; i++){
        x[i]=0.0;
    }
    //初期error
    std::vector<MY_FLOAT_TYPE> Ax_tmp(n);
    mat_vec_product(A, x, Ax_tmp, n);
    for(int ix=0;ix<n;ix++){
        error[ix]=b[ix]-Ax_tmp[ix];
//        direction[ix]=error[ix];
    }
    //初期directionの計算
    calc_LDLt_inv_r(L, d, error, direction, n);
    //初期エラーの大きさ
    MY_FLOAT_TYPE errerr_0=0.0;
    for(int ix=0;ix<n;ix++){
        errerr_0+=error[ix]*error[ix];
    }
    //メインの計算部分
    for(int i=0;i<max_itr;i++){
        mat_vec_product(A, direction, Adir, n);

        //解の修正係数coeff_0の計算
        std::vector<MY_FLOAT_TYPE> LDLt_inv_err(n);
        calc_LDLt_inv_r(L, d, error, LDLt_inv_err, n);
        //分子と分母を同時に計算
        MY_FLOAT_TYPE err_LDLt_inv_err_prev=0.0, diry=0.0;
        for(int ix=0;ix<n;ix++){
            err_LDLt_inv_err_prev+=error[ix]*LDLt_inv_err[ix];
            diry+=Adir[ix]*direction[ix];
        }
        MY_FLOAT_TYPE coeff_0=err_LDLt_inv_err_prev/diry;
        for(int ix=0;ix<n;ix++){
            //解の近似値を計算
            x[ix]+=coeff_0*direction[ix];
            //errorの計算
            error[ix]-=coeff_0*Adir[ix];
        }
        //ここで終了判定を行う
        //errorの大きさを計算
        MY_FLOAT_TYPE errerr_next=0.0;
        for(int ix=0;ix<n;ix++){
            errerr_next+=error[ix]*error[ix];
        }
        //終了条件
        if((errerr_next/errerr_0)<eps*eps){
/*
            std::ofstream writing_file;
            writing_file.open("iteration_num_ICCG_512_2d.dat", std::ios::app);
            writing_file << i << std::endl;
            writing_file.close();
*/
            return;
        }
/*
        if(errerr_next<eps*eps){
            return;
        }
*/

        //directionの修正係数coeff_0の計算
        calc_LDLt_inv_r(L, d, error, LDLt_inv_err, n);
        MY_FLOAT_TYPE err_LDLt_inv_err_next=0.0;
        for(int ix=0;ix<n;ix++){
            err_LDLt_inv_err_next+=error[ix]*LDLt_inv_err[ix];
        }
        MY_FLOAT_TYPE coeff_dir=err_LDLt_inv_err_next/err_LDLt_inv_err_prev;
        for(int ix=0;ix<n;ix++){
            direction[ix]=LDLt_inv_err[ix]+coeff_dir*direction[ix];
        }
    }
/*
    std::ofstream writing_file;
    writing_file.open("iteration_num_ICCG_512_2d.dat", std::ios::app);
    writing_file << max_itr << std::endl;
    writing_file.close();
*/
}

//gauss_seidel法sparseバージョン
void gauss_seidel(const sparse_matrix_with_diagonal_element &A, const std::vector<MY_FLOAT_TYPE> &b, std::vector<MY_FLOAT_TYPE> &x, int N, int max_itr){
    for(int ix=0;ix<N;ix++){
        x[ix]=0.0;
    }
    for(int k=0;k<max_itr;k++){
        int irow=0;
        int ielem=0;
        for(auto itr=A.cumulative_num_nonzero_element_in_row.begin()+1;itr!=A.cumulative_num_nonzero_element_in_row.end();itr++){
            MY_FLOAT_TYPE Ax_tmp=0.0;
            for(int i=0;i<(*(itr)-*(itr-1));i++){
                if(irow!=A.column_index[ielem]){
                    Ax_tmp+=A.element_value[ielem]*x[A.column_index[ielem]];
                }
                ielem++;
            }
            x[irow]=(b[irow]-Ax_tmp)/A.diagonal_element_value[irow];
            irow++;
        }
    }
}

}//namespace linear_algebra
