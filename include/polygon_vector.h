#pragma once

#include <vector>
#include <mutex>
#include "polygon_3D.h"
#include "define_float_type.h"

namespace smoke_simulation {
// マルチスレッドで使うためのpolygon_3Dのvector
// http://alpha-beta-gunma.hatenadiary.jp/entry/2018/01/29/181334 参照
class polygon_vector {
public:
//    std::mutex mtx;
    //data[i_thread][i_poly] は ithread 番目のスレッドの i_poly 番目のポリゴン
    std::vector<std::vector<polygon_3D>> data;
    const int num_threads;

    polygon_vector();
/*
    void push_back(polygon_3D x){
        std::lock_guard<std::mutex> lock(mtx);
        data.push_back(x);
    }
*/
};
}// namespace smoke_simulation
