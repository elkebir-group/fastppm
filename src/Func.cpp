//
// Created by Yuanyuan Qi on 6/26/24.
//

#include "Func.h"

//template<typename ctype>
int binary_search_less(const std::unordered_map<int, ctype> & sorted_list, int size_of_list, const ctype &value){
    int l_idx = 0, r_idx = size_of_list, mid_idx = (size_of_list)>>1;
    while (l_idx < r_idx - 1){//while size > 1
        if (sorted_list.at(mid_idx)<value){
            l_idx = mid_idx;
        } else {
            r_idx = mid_idx;
        }
        mid_idx = (l_idx+r_idx)>>1;
    }
    return mid_idx;
}

//template<typename ctype>
int binary_search_greater(const std::unordered_map<int, ctype> &sorted_list, int size_of_list, const ctype &value){
    int l_idx = 0, r_idx = size_of_list, mid_idx = (size_of_list)>>1;
    while (l_idx < r_idx - 1){//while size > 1
        if (sorted_list.at(mid_idx)>value){
            l_idx = mid_idx;
        } else {
            r_idx = mid_idx;
        }
        mid_idx = (l_idx+r_idx)>>1;
    }
    return mid_idx;
}

real log_eps(real x, real eps, int s_n) {
    if (x < eps) {
        real return_val = log(eps);
        real iter_val = (eps - x) / eps;
        for (int i=1; i < s_n; i++) {
            return_val -= pow(iter_val , i) / i;
        }
        return return_val;
    }
    return log(x);
}
