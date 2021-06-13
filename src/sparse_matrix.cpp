#include "sparse_matrix.h"

#include "define_float_type.h"

namespace linear_algebra{

sparse_matrix::sparse_matrix(int nr, int nc){
    num_element=0;
    num_row=nr;
    num_column=nc;
    input_row_num=0;
    //input_element関数での帳尻合わせ
    cumulative_num_nonzero_element_in_row.push_back(0);
    cumulative_num_nonzero_element_in_row.push_back(0);
}

//row行目、column列目の値elementを入力する関数
void sparse_matrix::input_element(int row, int column, MY_FLOAT_TYPE element){
    if(row<input_row_num){
        std::cerr<<"入れる行が戻っている  row: "<<row<<",  column: "<<column<<std::endl;
    }
    else{
        if(row==input_row_num){
            element_value.push_back(element);
            column_index.push_back(column);
            cumulative_num_nonzero_element_in_row[input_row_num+1]++;
        }
        else{
            //空の行を埋める
            while(row!=input_row_num){
                cumulative_num_nonzero_element_in_row.push_back(cumulative_num_nonzero_element_in_row[input_row_num+1]);
                input_row_num++;
            }
            element_value.push_back(element);
            column_index.push_back(column);
            cumulative_num_nonzero_element_in_row[input_row_num+1]++;
        }
    }
}
//dense形式で行列を表示する関数
void sparse_matrix::print_matrix(){
    int ielem=0;
    for(int ix=0;ix<num_row;ix++){
        int icol=0;
        for(int iy=0;iy<num_column;iy++){
            if(icol<(cumulative_num_nonzero_element_in_row[ix+1]-cumulative_num_nonzero_element_in_row[ix])){
                if(iy==column_index[ielem]){
                    std::cout<<element_value[ielem]<<" ";
                    ielem++;
                    icol++;
                }
                else{
                    std::cout<<0<<" ";
                }
            }
            else{
                std::cout<<0<<" ";
            }
        }
        std::cout<<std::endl;
    }
}

void sparse_matrix::print_information(){
    std::cout<<"element_value: "<<std::endl;
    for(auto itr=element_value.begin();itr!=element_value.end();itr++){
        std::cout<<*itr<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"column_index: "<<std::endl;
    for(auto itr=column_index.begin();itr!=column_index.end();itr++){
        std::cout<<*itr<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"cumulative_num_nonzero_element_in_row: "<<std::endl;
    for(auto itr=cumulative_num_nonzero_element_in_row.begin();itr!=cumulative_num_nonzero_element_in_row.end();itr++){
        std::cout<<*itr<<" ";
    }
    std::cout<<std::endl;
    std::cout<<std::endl;
}


sparse_matrix_with_diagonal_element::sparse_matrix_with_diagonal_element(int nr, int nc): sparse_matrix(nr, nc){
    int n_max=std::max(nr,nc);
    diagonal_element_value.resize(n_max);
    for(int i=0;i<n_max;i++){
        diagonal_element_value[i]=0.0;
    }
}

void sparse_matrix_with_diagonal_element::input_element(int row, int column, MY_FLOAT_TYPE element){
    if(row<input_row_num){
        std::cerr<<"入れる行が戻っている  row: "<<row<<",  column: "<<column<<std::endl;
    }
    else{
        if(row==input_row_num){
            element_value.push_back(element);
            column_index.push_back(column);
            cumulative_num_nonzero_element_in_row[input_row_num+1]++;
            if(row==column){
                diagonal_element_value[row]=element;
            }
        }
        else{
            while(row!=input_row_num){
                cumulative_num_nonzero_element_in_row.push_back(cumulative_num_nonzero_element_in_row[input_row_num+1]);
                input_row_num++;
            }
            element_value.push_back(element);
            column_index.push_back(column);
            cumulative_num_nonzero_element_in_row[input_row_num+1]++;
        }
    }
}

//sparse_matrixとvector<MY_FLOAT_TYPE>の積を計算する。結果は引数で渡されたxに格納される
void mat_vec_product(const sparse_matrix A,const std::vector<MY_FLOAT_TYPE> b, std::vector<MY_FLOAT_TYPE>& x, int n){
    for(int i=0;i<n;i++){
        x[i]=0.0;
    }
    int irow=0;
    int ielem=0;
    //全ての行を走査するループ
    for(auto itr=A.cumulative_num_nonzero_element_in_row.begin()+1;itr!=A.cumulative_num_nonzero_element_in_row.end();itr++){
        //各行における非ゼロな値を持つ列を舐めるループ
        for(int i=0;i<(*(itr)-*(itr-1));i++){
            x[irow]+=A.element_value[ielem]*b[A.column_index[ielem]];
            ielem++;
        }
        irow++;
    }
}

}//namespace linear_algebra
