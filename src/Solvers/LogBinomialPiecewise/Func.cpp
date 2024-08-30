//
// Created by Yuanyuan Qi on 6/26/24.
//

#include "Func.h"

namespace LogBinomialPiecewiseLinearSolver {

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

real func_llh(real para_var, real para_ref, real x) {
    return -para_var * log_eps(x) - para_ref * log_eps(1 - x);
}

};
