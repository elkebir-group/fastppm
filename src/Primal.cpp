//
// Created by Yuanyuan Qi on 8/1/24.
//

#include "Primal.h"

PWL_close::PWL_close(int max_k) :
        k (0),
        x(max_k + 1),
        y(max_k + 1),
        slope(max_k) {
}

void PWL_close::update(real begin, real end, int _k, const std::pair<real,real> &func_para,
                       std::function<real(real,real,real)> & lossFunction) {
    k = _k;
    for (int i = 0; i <= k; i++) {
        x[i] = begin + (end - begin) * i / k;
        y[i] = lossFunction(func_para.first, func_para.second, x[i]);
    }
    for (int i = 0; i < k; i++) {
        slope[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    }
}

real PWL_close::backtrace(real d, real * debug){
    int idx = std::upper_bound(slope.begin(),slope.begin()+k,d)-slope.begin();
    if (idx>0 && (d-slope[idx-1]) < Compare_eps) idx--;
    if (debug!=NULL){
        *debug = y[idx] - x[idx]*d;
    }
    return x[idx];
}

#ifdef _DEBUG

real PWL_close::operator()(real xx) const {
    int idx = std::lower_bound(x.begin(),x.begin()+k,xx)-x.begin();
    if (idx <=0) return y[0];
    return y[idx]+(xx-x[idx])*slope[idx-1];
}


bool PWL_close::self_check() const {
    for (int i = 0; i < k; i++) {
        if (x[i+1] < x[i]){
            return false;
        }
        if (abs((y[i+1]-y[i])/(x[i+1]-x[i]) -slope[i]) > Compare_eps){
            return false;
        }
    }
    return true;
}

#endif