//
// Created by Yuanyuan Qi on 6/24/24.
//

#ifndef EFFICIENTLLHESTIMATOR_PWL_H
#define EFFICIENTLLHESTIMATOR_PWL_H

#include <functional>
#include <vector>

#include "Func.h"

class PWL { //convex, obj, min
public:
    PWL(int n_intervals);

    // objective_func // convex function // must not be increasing
    int k;
    //linear pw objective_func //convex function
    std::vector<real> prime_x; // k+1
    std::vector<real> prime_y; // k+1
    std::vector<real> prime_slope; // k //increasing

    //objective with dual variables :min over original variable, max over dual //decreasing and concave function
    //d (delta) = a_i - a_pi (pi: parent of i)// a_proot = beta
//    std::vector<real> & dual_d_x; //k //identical to prime_slope
    std::vector<real> dual_d_y;  //k
//    std::vector<real> & neg_dual_d_slope;//k+1 //dual_d_slope is identical to -prime_x which is decreasing

    void update(real begin, real end, const Func & func);

    real bt_f(real d_dual);
};

#endif //EFFICIENTLLHESTIMATOR_PWL_H
