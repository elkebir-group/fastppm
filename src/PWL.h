//
// Created by Yuanyuan Qi on 6/24/24.
//

#ifndef EFFICIENTLLHESTIMATOR_PWL_H
#define EFFICIENTLLHESTIMATOR_PWL_H

#include <functional>
#include <vector>

#include "Func.h"

class PWL_close;
typedef PWL_close Primal;

class PWL_open;
typedef PWL_open Dual;

class PWL_Ropen{// x bounded on left side // left is always 0. i.e. x[0]=0 // state function
public:
    PWL_Ropen(int max_k);
    int k; // n of intervals;
    std::vector<real> x;
    std::vector<real> y;
    std::vector<real> slope;

    void base_case(const Dual & dual);
    void sum(const std::vector<PWL_Ropen *> & to_sum, std::vector<std::pair<real, real> > & helper);
    void optimize_with(const PWL_Ropen & sum_func, const Dual & dual);
    real backtrace(const Dual &dual, real val, std::vector<real> & helper_x, std::vector<real> & helper_y) const ;

#ifdef _DEBUG
    real operator()(real x) const;
    bool self_check() const;
#endif
};

//
//class PWL { //convex, obj, min
//public:
//    PWL(int n_intervals);
//
//    // objective_func // convex function // must not be increasing
//    int k;
//    //linear pw objective_func //convex function
//    std::vector<real> prime_x; // k+1
//    std::vector<real> prime_y; // k+1
//    std::vector<real> prime_slope; // k //increasing
//
//    //objective with dual variables :min over original variable, max over dual //decreasing and concave function
//    //d (delta) = a_i - a_pi (pi: parent of i)// a_proot = beta
////    std::vector<real> & dual_d_x; //k //identical to prime_slope
//    std::vector<real> dual_d_y;  //k
////    std::vector<real> & neg_dual_d_slope;//k+1 //dual_d_slope is identical to -prime_x which is decreasing
//
//    void update(real begin, real end, const Func & func);
//
//    real bt_f(real d_dual);
//};

#endif //EFFICIENTLLHESTIMATOR_PWL_H
