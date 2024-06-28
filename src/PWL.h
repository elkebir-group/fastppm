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
};

class DPSTATE{ //concave and increasing, max, w r t to a single var, defined on [0,inf)
public:
    //result
    int effective_k;
    std::vector<real> x; //k
    std::vector<real> y;
    std::vector<real> slope; //k+1 //decreasing

    //sum //concave and increasing
    int cnt, nbrk;
    std::vector<std::pair<real, real> > brk_points;
    std::vector<real> sum_x;
    std::vector<real> sum_y;
    std::vector<real> sum_slope; // decreasing

    std::vector<real>::iterator it;//tmp_var

    DPSTATE();

    void init(int n_children, int k);

    void zero(PWL &a);

    void update(std::vector<DPSTATE*> & to_sum, PWL & a);
};

#endif //EFFICIENTLLHESTIMATOR_PWL_H
