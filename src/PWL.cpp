//
// Created by Yuanyuan Qi on 6/24/24.
//

#include "PWL.h"
#include "Func.h"

PWL::PWL(int n_intervals) :
        k(n_intervals),
        prime_x(n_intervals + 1),
        prime_y(n_intervals + 1),
        prime_slope(n_intervals),
//        dual_d_x(prime_slope),
        dual_d_y(n_intervals)
//        neg_dual_d_slope(prime_x)
{
}

void PWL::update(real begin, real end, const Func & func) {
    for (int i = 0; i <= k; i++) {
        prime_x[i] = begin + (end - begin) * i / k;
        prime_y[i] = -func(prime_x[i]);
    }
    for (int i = 0; i < k; i++) {
        prime_slope[i] = (prime_y[i + 1] - prime_y[i]) / (prime_x[i + 1] - prime_x[i]);
        //note that dual_d_x = prime_slope;
        dual_d_y[i] = prime_y[i] - prime_slope[i] * prime_x[i]; //actually intercept
    }
    //dual_d_slope=-prime_x
}

#define BT_eps 1e-7
real PWL::bt_f(real d_dual) {
    int idx = std::upper_bound(prime_slope.begin(),prime_slope.end(),d_dual)-prime_slope.begin();
    if (idx>0 && (d_dual-prime_slope[idx-1])<BT_eps) idx--;
    return prime_x[idx];
}
#undef BT_eps

