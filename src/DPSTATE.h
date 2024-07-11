//
// Created by Yuanyuan Qi on 7/11/24.
//

#ifndef EFFICIENTLLHESTIMATOR_DPSTATE_H
#define EFFICIENTLLHESTIMATOR_DPSTATE_H

#include <vector>
#include "Func.h"

class PWL;

class DPSTATE{ //concave and increasing, max, w r t to a single var, defined on [0,inf)
public:
    //result
    int effective_k;
    std::vector<real> x; //k+1 (0 included)
    std::vector<real> y;
    std::vector<real> slope; //k+1 //decreasing
    std::vector<real> prime_x_BT;
    std::vector<real> dual_d_a_BT;

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

    real traceback_x;
    real backtrace(real val);
};

#endif //EFFICIENTLLHESTIMATOR_DPSTATE_H
