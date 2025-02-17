//
// Created by Yuanyuan Qi on 6/25/24.
//

#ifndef EFFICIENTLLHESTIMATOR_FUNC_H
#define EFFICIENTLLHESTIMATOR_FUNC_H
#include <cmath>
#include <vector>
#include <unordered_map>
#include <iostream>

namespace PiecewiseLinearSolver {

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

inline real FuncBinomLlh(real para_var, real para_ref, real x) {
    if (para_var == 0)
        return -para_ref * log_eps(1 - x);
    else if (para_ref == 0)
        return -para_var * log_eps(x);
    else
        return -para_var * log_eps(x) - para_ref * log_eps(1 - x);
}

inline real FuncBetaBinomLlh(real para_precision, real para_var, real para_ref, real x) {
//     \log\Gamma(f_i s)-\log\Gamma(a_i + f_i s) + \log\Gamma(s-f_i s) - \log\Gamma(d_i-a_i+s -f_i s)
//    real ans = std::lgamma(x * para_precision) - std::lgamma(para_var + x * para_precision)
//               + std::lgamma(para_precision - x * para_precision) - std::lgamma(para_ref + para_precision - x * para_precision);
//    std::cout << para_precision << " " << para_var << " " << para_ref << " " << x << " " << ans << std::endl;
//
//    if (std::isnan(ans))
//    {
//        std::cerr << "break" << std::endl;
//    }
//
//    return ans;
    return std::lgamma(x * para_precision) - std::lgamma(para_var + x * para_precision)
        + std::lgamma(para_precision - x * para_precision) - std::lgamma(para_ref + para_precision - x * para_precision);
}

class func_base {
public:
    real var, ref;
    virtual real operator()(real x) = 0;
};

class func_binom_llh: public func_base {
public:
    inline real operator()(real x) override {
        return FuncBinomLlh(var, ref, x);
    }
};

class func_beta_binom_llh: public func_base {
public:
    static real precision;
    inline real operator()(real x) override {
        return FuncBetaBinomLlh(precision, var, ref, x);
    }
};

};

#endif //EFFICIENTLLHESTIMATOR_FUNC_H
