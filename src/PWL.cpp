//
// Created by Yuanyuan Qi on 6/24/24.
//

#include "PWL.h"
#include "Func.h"
#include <queue>
#include <cassert>

typedef std::pair<PWL_Ropen*, int> priority_queue_helper;
inline bool operator < (const priority_queue_helper & a, const priority_queue_helper & b){
    return a.first->x[a.second] > b.first->x[b.second];
}

PWL_close::PWL_close(int max_k) :
        k (0),
        x(max_k + 1),
        y(max_k + 1),
        slope(max_k) {
}

void PWL_close::update(real begin, real end, int _k, const Func &func) {
    k = _k;
    for (int i = 0; i <= k; i++) {
        x[i] = begin + (end - begin) * i / k;
        y[i] = func(x[i]);
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

PWL_open::PWL_open(int max_k) :
        k (0),
        x(max_k - 1),
        y(max_k - 1),
        slope(max_k) {
}

void PWL_open::dual_from_primal(const PWL_close &primal) {
    k = primal.k + 1;
    for (int i = 0; i < k - 1; i++) {
        x[i] = primal.slope[i];
        y[i] = primal.y[i] - primal.slope[i] * primal.x[i];
    }
    for (int i = 0; i < k; i++) {
        slope[i] = -primal.x[i];
    }
}

real PWL_close::operator()(real xx) const {
    int idx = std::lower_bound(x.begin(),x.begin()+k,xx)-x.begin();
    if (idx <=0) return y[0];
    return y[idx]+(xx-x[idx])*slope[idx-1];
}

real PWL_open::operator()(real xx) const {
    int idx = std::lower_bound(x.begin(),x.begin()+k-1,xx)-x.begin();
    return y[idx]+(xx-x[idx])*slope[idx];
}

real PWL_Ropen::operator()(real xx) const {
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

bool PWL_open::self_check() const {
    for (int i = 0; i < k - 2; i++) {
        if (x[i+1] <= x[i]){
            return false;
        }
        if (abs((y[i+1]-y[i])/(x[i+1]-x[i]) - slope[i+1]) > Compare_eps){
            return false;
        }
    }
    return true;
}

bool PWL_Ropen::self_check() const {
    for (int i = 0; i < k - 1; i++) {
        if (x[i+1] <= x[i]){
            return false;
        }
        if (abs((y[i+1]-y[i])/(x[i+1]-x[i]) - slope[i]) > assert_eps){
            printf("%.12lf\n",(y[i+1]-y[i])/(x[i+1]-x[i]) - slope[i]);
            printf("%d %.12lf %.12lf %.12lf %.12lf %.12lf R\n",i,y[i+1],y[i],x[i+1],x[i], slope[i]);
            return false;
        }
    }
    return true;
}


PWL_Ropen::PWL_Ropen(int max_k) :
        k (0),
        x(max_k),
        y(max_k),
        slope(max_k){
}

void PWL_Ropen::base_case(const PWL_open & dual) {
    k = std::upper_bound(dual.x.begin(),dual.x.begin()+dual.k - 1,0.) - dual.x.begin();
    if (k>=dual.k-1) y[0] = dual.y[k-1] - dual.x[k-1] * dual.slope[k];
    else y[0] = dual.y[k] - dual.x[k] * dual.slope[k];
    x[0] = 0;
    for (int i = 0; i < k; i++){
        x[i+1] = -dual.x[k - i - 1];
        y[i+1] = dual.y [k - i - 1];
    }
    for (int i = 0; i <= k; i++){
        slope[i] = -dual.slope[k - i];
    }
    k++;
}

void PWL_Ropen::sum(const std::vector<PWL_Ropen *> & to_sum, std::vector<std::pair<real, real> > & helper) {
    int brk_idx = 0;
    x[0] = 0;
    y[0] = 0;
    slope[0] = 0;
    std::priority_queue<priority_queue_helper> _priority_queue;

    for (auto it: to_sum) {
        _priority_queue.emplace(priority_queue_helper(it, 1));
        slope[0] += it->slope[0];
        y[0] += it->y[0];
    }

    priority_queue_helper _top;
    while (!_priority_queue.empty()) {
        _top = _priority_queue.top();
        _priority_queue.pop();
        if (_top.second >= _top.first->k) continue;
        helper[brk_idx].first = _top.first->x[_top.second];
        helper[brk_idx].second = _top.first->slope[_top.second] - _top.first->slope[_top.second - 1];
        brk_idx++;
        _top.second++;
        _priority_queue.push(_top);
    }

    k = 1;
    for (int i = 0; i < brk_idx; k++) {
        slope[k] = slope[k - 1];
        x[k] = helper[i].first;
        while (i < brk_idx && helper[i].first <= x[k]) {
            slope[k] += helper[i].second;
            i++;
        }
    }

    for (int i = 1; i < k; i++) {
        y[i] = y[i - 1] + (x[i] - x[i - 1]) * slope[i - 1];
    }
}

void PWL_Ropen::optimize_with(const PWL_Ropen &sum_func, const PWL_open &dual) {
    if (dual.slope[dual.k-1] + sum_func.slope[sum_func.k-1] > Compare_eps){ // note 1;
        printf("ERROR, oh no!\n");
    }
    if (dual.slope[0] + sum_func.slope[0] <= Compare_eps ){ //note 2; decreasing always, easy
        base_case(dual);
        for (int i = 0; i < k; i++) {
            y[i] += sum_func.y[0];
        }
        return;
    }

    //increase and then decrease over a_i
    k = 0;
    int dual_idx = 0, sum_idx = sum_func.k-1;
    if (dual.slope[dual_idx] + sum_func.slope[sum_idx] <= 0){ //0 cnt
        while (dual.slope[dual_idx] + sum_func.slope[sum_idx-1] <= 0) sum_idx --; // will stop because of note 2
    } else {
        while (dual.slope[dual_idx+1] + sum_func.slope[sum_idx] > 0) dual_idx++; // will stop because of note 1
    }

    //slope last
    if (dual.slope[dual_idx] + sum_func.slope[sum_idx] <= 0){
        slope[0] = -dual.slope[dual_idx];
    } else {
        slope[0] = sum_func.slope[sum_idx];
    }

    if (sum_func.x[sum_idx] - dual.x[dual_idx] <= 0){
        x[0] = sum_func.x[sum_idx] - dual.x[dual_idx];
        y[0] = dual.y[dual_idx] + sum_func.y[sum_idx];
        y[0] = y[0] + (0 - x[0]) * slope[0];
        x[0] = 0;
    }
    else {
        while (dual_idx < dual.k-1 && sum_func.x[sum_idx] - dual.x[dual_idx] > 0) { //sumidx >=0
            x[k] = sum_func.x[sum_idx] - dual.x[dual_idx];
            y[k] = sum_func.y[sum_idx] + dual.y[dual_idx];
            //real
            k++;
            if (sum_idx >= 1 && dual.slope[dual_idx + 1] + sum_func.slope[sum_idx - 1] <= 0) {
                sum_idx--;
                slope[k] = sum_func.slope[sum_idx];
            } else {
                dual_idx++;
                slope[k] = -dual.slope[dual_idx];
            }
        }

        x[k] = 0;
        y[k] = y[k - 1] - x[k - 1] * slope[k];
    }
    k++;
    std::reverse(x.begin(), x.begin() + k);
    std::reverse(y.begin(), y.begin() + k);
    std::reverse(slope.begin(), slope.begin() + k);
    real y_re;
    for (int i = 0; i < k; i++){
        sum_func.backtrace(dual,x[i], &y_re, false);
        if (abs(y_re-y[i]) > Compare_eps){
            sum_func.backtrace(dual,x[i], &y_re,true);
            assert(abs(y_re-y[i]) < Compare_eps);
        }
    }
}


real PWL_Ropen::backtrace(const PWL_open &dual, real val, real * debug, bool flag) const {
    std::vector<real> result_x (k + dual.k + 2), result_y(k + dual.k + 2);
    int sum_idx=0, pwl_idx=0, tot_index=0;
    real _x;
    while (sum_idx < k && pwl_idx < dual.k - 1) {
        _x = dual.x[pwl_idx] + val;
        if (_x < 0) {
            pwl_idx++;
            continue;
        }
        if (x[sum_idx] == _x) {
            result_x[tot_index] = _x;
            result_y[tot_index] = y[sum_idx] + dual.y[pwl_idx];
            pwl_idx++;
            sum_idx++;
        } else if (x[sum_idx] < _x) {
            result_x[tot_index] = x[sum_idx];
            result_y[tot_index] = y[sum_idx] + dual.y[pwl_idx] + (x[sum_idx] - _x) * dual.slope[pwl_idx];
            sum_idx++;
        } else {
            result_x[tot_index] = _x;
            result_y[tot_index] = dual.y[pwl_idx] + y[sum_idx] + (_x - x[sum_idx]) * slope[sum_idx - 1];
            pwl_idx++;
        }
        tot_index++;
    }
    while (sum_idx < k) {
        result_x[tot_index] = x[sum_idx];
        result_y[tot_index] =
                y[sum_idx] + dual.y[dual.k - 2] + (x[sum_idx] - dual.x[dual.k - 2] - val) * dual.slope[dual.k - 1];
        sum_idx++;
        tot_index++;
    }
    while (pwl_idx < dual.k-1) {
        _x = dual.x[pwl_idx] + val;
        result_x[tot_index] = _x;
        result_y[tot_index] = dual.y[pwl_idx] + y[k - 1] + (_x - x[k - 1]) * slope[k - 1];
        pwl_idx++;
        tot_index++;
    }

//    printf("-------------------------------\n");
    real xx=0, yy=-1e300;

    for (int i = 0; i < tot_index; i++){
        if (result_y[i] - yy > 0){
            xx = result_x[i];
            yy = result_y[i];
        }
        if (flag)printf("%lf %lf\n",xx,yy);
    }
    if (debug)*debug = yy;
    return xx;
}
