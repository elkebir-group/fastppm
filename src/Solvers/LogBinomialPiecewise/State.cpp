//
// Created by Yuanyuan Qi on 6/24/24.
//


#include <queue>
#include <algorithm>
#ifdef _DEBUG
#include <cassert>
#endif

#include "State.h"
#include "Func.h"
#include "Primal.h"
#include "Dual.h"

namespace LogBinomialPiecewiseLinearSolver {
std::unordered_map<int, std::pair<real,real> > State::helper;
std::unordered_map<int, real> State::helper_x, State::helper_y;

typedef std::pair<const State*, int> priority_queue_helper;
inline bool operator < (const priority_queue_helper & a, const priority_queue_helper & b) {
    return a.first->x.at(a.second) > b.first->x.at(b.second);
}

void State::base_case(const step_base * dual_ptr) {
    const Dual & dual (*(Dual *) dual_ptr);
//    k = binary_search_less(dual.x, dual.k - 1, 0.);
    k = std::upper_bound(dual.x.begin(),dual.x.begin()+dual.k - 1,0.) - dual.x.begin();
    if (k>=dual.k-1) y[0] = dual.y[k-1] - dual.x[k-1] * dual.slope[k];
    else y[0] = dual.y[k] - dual.x[k] * dual.slope[k];
    x[0] = 0;
    for (int i = 0; i < k; i++){
        x[i+1] = -dual.x[k - i - 1];
        y[i+1] = dual.y[k - i - 1];
    }
    for (int i = 0; i <= k; i++){
        slope[i] = -dual.slope[k - i];
    }
    k++;
}

void State::sum(const std::vector<const state_base *> & to_sum) {
    int brk_idx = 0;
    x[0] = 0;
    y[0] = 0;
    slope[0] = 0;
    std::priority_queue<priority_queue_helper> _priority_queue;

    for (auto it: to_sum) {
        if (((const State *)it)->k >= 2)
            _priority_queue.emplace(priority_queue_helper((const State *)it, 1));
        slope[0] += ((const State *)it)->slope.at(0);
        y[0] += ((const State *)it)->y.at(0);
    }

    priority_queue_helper _top;
    while (!_priority_queue.empty()) {
        _top = _priority_queue.top();
        _priority_queue.pop();
        helper[brk_idx].first = _top.first->x.at(_top.second);
        helper[brk_idx].second = _top.first->slope.at(_top.second) - _top.first->slope.at(_top.second - 1);
        brk_idx++;
        _top.second++;
        if (_top.second >= _top.first->k) continue;
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

void State::optimize_with(const sum_base *sum_func_ptr, const step_base * dual_ptr) {
    const State &sum_func(*(State *)sum_func_ptr);
    const Dual &dual(*(Dual *) dual_ptr);
#ifdef _DEBUG
    assert (dual.slope.at(dual.k-1) + sum_func.slope.at(sum_func.k-1) < assert_eps); // note 1;
#endif
    if (dual.slope[0] + sum_func.slope.at(0) <= Compare_eps ){ //note 2; decreasing always, easy
        base_case(dual_ptr);
        for (int i = 0; i < k; i++) {
            y[i] += sum_func.y.at(0);
        }
        return;
    }

    //increase and then decrease over a_i
    k = 0;
    int dual_idx = 0, sum_idx = sum_func.k-1;
    if (dual.slope[dual_idx] + sum_func.slope.at(sum_idx) <= 0){ //0 cnt
        while (dual.slope[dual_idx] + sum_func.slope.at(sum_idx-1) <= 0) { // will stop because of note 2
            sum_idx--;
        }
    } else {
        while (dual.slope[dual_idx + 1] + sum_func.slope.at(sum_idx) > 0) {// will stop because of note 1
            //printf("%.6lf\n",dual.slope.at(dual_idx + 1) + sum_func.slope.at(sum_idx));
            dual_idx++;
        }
    }

    //slope last
    if (dual.slope[dual_idx] + sum_func.slope.at(sum_idx) <= 0){
        slope[0] = -dual.slope[dual_idx];
    } else {
        slope[0] = sum_func.slope.at(sum_idx);
    }

    if (sum_func.x.at(sum_idx) - dual.x[dual_idx] <= 0){
        x[0] = sum_func.x.at(sum_idx) - dual.x[dual_idx];
        y[0] = dual.y[dual_idx] + sum_func.y.at(sum_idx);
        y[0] = y[0] + (0 - x[0]) * slope[0];
        x[0] = 0;
    }
    else {
        while (dual_idx < dual.k-1 && sum_func.x.at(sum_idx) - dual.x[dual_idx] > 0) { //sumidx >=0
            x[k] = sum_func.x.at(sum_idx) - dual.x[dual_idx];
            y[k] = sum_func.y.at(sum_idx) + dual.y[dual_idx];
            //real
            k++;
            if (sum_idx >= 1 && dual.slope[dual_idx + 1] + sum_func.slope.at(sum_idx - 1) <= 0) {
                sum_idx--;
                slope[k] = sum_func.slope.at(sum_idx);
            } else {
                dual_idx++;
                slope[k] = -dual.slope[dual_idx];
            }
        }

        x[k] = 0;
        y[k] = y[k - 1] - x[k - 1] * slope[k];
    }
    k++;
    for (int i = 0; i < (k>>1); i++){
        std::swap(x[i],x[k-1-i]);
        std::swap(y[i],y[k-1-i]);
        std::swap(slope[i], slope[k-1-i]);
    }
}

real State::backtrace(const step_base *dual_ptr, real val) const {
    const Dual &dual(*(Dual*) dual_ptr);
    int sum_idx = 0, pwl_idx = 0, tot_index = 0;
    real _x;
    while (sum_idx < k && pwl_idx < dual.k - 1) {
        _x = dual.x[pwl_idx] + val;
        if (_x < 0) {
            pwl_idx++;
            continue;
        }
        if (x.at(sum_idx) == _x) {
            helper_x[tot_index] = _x;
            helper_y[tot_index] = y.at(sum_idx) + dual.y[pwl_idx];
            pwl_idx++;
            sum_idx++;
        } else if (x.at(sum_idx) < _x) {
            helper_x[tot_index] = x.at(sum_idx);
            helper_y[tot_index] = y.at(sum_idx) + dual.y[pwl_idx] + (x.at(sum_idx) - _x) * dual.slope[pwl_idx];
            sum_idx++;
        } else {
            helper_x[tot_index] = _x;
            helper_y[tot_index] = dual.y[pwl_idx] + y.at(sum_idx) + (_x - x.at(sum_idx)) * slope.at(sum_idx - 1);
            pwl_idx++;
        }
        tot_index++;
    }
    while (sum_idx < k) {
        helper_x[tot_index] = x.at(sum_idx);
        helper_y[tot_index] =
                y.at(sum_idx) + dual.y[dual.k - 2] + (x.at(sum_idx) - dual.x[dual.k - 2] - val) * dual.slope[dual.k - 1];
        sum_idx++;
        tot_index++;
    }
    while (pwl_idx < dual.k - 1) {
        _x = dual.x[pwl_idx] + val;
        helper_x[tot_index] = _x;
        helper_y[tot_index] = dual.y[pwl_idx] + y.at(k - 1) + (_x - x.at(k - 1)) * slope.at(k - 1);
        pwl_idx++;
        tot_index++;
    }

    real xx = 0, yy = -1e300;

    for (int i = 0; i < tot_index; i++) {
        if (helper_y[i] - yy > 0) {
            xx = helper_x[i];
            yy = helper_y[i];
        }
    }
    return xx;
}

#ifdef _DEBUG

real State::operator()(real xx) const {
    int idx = std::lower_bound(x.begin(),x.begin()+k,xx)-x.begin();
    if (idx <=0) return y[0];
    return y[idx-1]+(xx-x[idx-1])*slope[idx-1];
}

bool State::self_check() const {
    for (int i = 0; i < k - 1; i++) {
        if (x[i+1] <= x[i]){
            return false;
        }
        if (abs((y[i+1]-y[i])/(x[i+1]-x[i]) - slope[i]) > assert_eps){
            return false;
        }
    }
    return true;
}

#endif

};

