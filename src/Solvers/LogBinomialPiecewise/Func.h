//
// Created by Yuanyuan Qi on 6/25/24.
//

#ifndef EFFICIENTLLHESTIMATOR_FUNC_H
#define EFFICIENTLLHESTIMATOR_FUNC_H
#include <cmath>
#include <vector>
#include <unordered_map>

namespace LogBinomialPiecewiseLinearSolver {

typedef double real;

#ifdef _DEBUG
#include <cassert>
const real assert_eps = 1.5e-6;
#endif

const real Compare_eps = 1e-7;

//template<typename ctype>
typedef real ctype;
int binary_search_less(const std::unordered_map<int, ctype> &sorted_list, int size_of_list, const ctype &value);

//template<typename ctype>
int binary_search_greater(const std::unordered_map<int, ctype> &sorted_list, int size_of_list, const ctype &value);

real log_eps(real x,real eps=1e-6,int s_n=3);

inline real Funcllh(real para_var, real para_ref, real x){
    return -para_var * log_eps(x) - para_ref * log_eps(1 - x);
}

class func_base{
public:
    virtual real operator()(real x) = 0;
};

class func_llh: public func_base{
public:
    real var,ref;
    inline real operator()(real x) override{
        return Funcllh(var, ref, x);
    }
};
};

#endif //EFFICIENTLLHESTIMATOR_FUNC_H
