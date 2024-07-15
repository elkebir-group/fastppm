//
// Created by Yuanyuan Qi on 7/11/24.
//

#include "DPSTATE.h"
#include "PWL.h"

#define dual_d_x prime_slope
#define neg_dual_d_slope prime_x

DPSTATE::DPSTATE()
{

}

void DPSTATE::init(int n_descendant, int k) {
    x.resize(k * (n_descendant + 1) + 2);
    y.resize(k * (n_descendant + 1) + 2);
    slope.resize(k * (n_descendant + 1) + 2);
    brk_points.resize(n_descendant * k + 1);
    prime_x_BT.resize(k * (n_descendant + 1) + 2);
    dual_d_a_BT.resize(k * (n_descendant + 1) + 2);
    sum_x.resize(n_descendant * k + 1);
    sum_y.resize(n_descendant * k + 1);
    sum_slope.resize(n_descendant * k + 2);
}

void DPSTATE::update(std::vector<DPSTATE*> & to_sum, PWL &a, bool debug) {
    if (to_sum.empty()) {//leaf a_i = 0, h(-d)
        zero(a);
        return;
    }

    // improvable with heap: log(nk) -> log(n) not really necessary
    nbrk = 0;
    sum_x[0] = 0;
    sum_y[0] = 0;
    sum_slope[0] = 0;
    for (auto it: to_sum) {
        for (int i = 1; i <= it->effective_k; i++) {
            brk_points[nbrk].first = it->x[i];
            brk_points[nbrk].second = it->slope[i] - it->slope[i-1];
            nbrk++;
        }
        sum_slope[0] += it->slope[0];
        sum_y[0] += it->y[0];
    }
    std::sort(brk_points.begin(), brk_points.begin() + nbrk);

    cnt = 0;
    for (int i = 0; i < nbrk; ) {
        cnt ++;
        sum_slope[cnt] = sum_slope[cnt-1];
        sum_x[cnt] = brk_points[i].first;
        while(i <= nbrk && brk_points[i].first <= sum_x[cnt]) {
            sum_slope[cnt] += brk_points[i].second;
            i++;
        }
    }

    for (int i = 1; i <= cnt; i++) {
        sum_y[i] = sum_y[i - 1] + (sum_x[i] - sum_x[i - 1]) * sum_slope[i-1];
    }

    //add h(d) and max over a_i, never strictly increasing => a neg d slop[-1] >= sum_slope[-1]
    if (a.neg_dual_d_slope[a.k] < sum_slope[cnt] || debug){
        printf("We have a problem!! %lf %lf\n",a.neg_dual_d_slope[a.k],sum_slope[cnt]);
    }
    if (a.neg_dual_d_slope[0] >= sum_slope[0]) { //always decreasing // a_i = 0 <=> lb of [fp] >= s ub of [fc]
        zero(a);
        y[0] += sum_y[0];
        return;
    }
    //increase and then decrease over a_i
    effective_k = 0;
    int aidx = 0, sumidx = cnt;
    bool break_flag = false;
    //start from +inf
    while (true) {
        if (a.neg_dual_d_slope[aidx] >= sum_slope[sumidx]) {
            while (sumidx >= 1 && sum_slope[sumidx - 1] <= a.neg_dual_d_slope[aidx]) {
                sumidx--;
            }
            while (aidx <= a.k && (sumidx <= 0 || a.neg_dual_d_slope[aidx] < sum_slope[sumidx - 1])) {
                //a_i = sum_x[sumidx], ([idx-1]<)a_i - a_pi < d_x[aidx]
                if (aidx >= a.k  || sum_x[sumidx] <= a.dual_d_x[aidx]) {
                    dual_d_a_BT[effective_k] = sum_x[sumidx];
                    slope[effective_k] = a.neg_dual_d_slope[aidx];
                    break_flag = true;
                    break;
                }
                dual_d_a_BT[effective_k] = sum_x[sumidx];
                x[effective_k] = sum_x[sumidx] - a.dual_d_x[aidx];
                y[effective_k] = sum_y[sumidx] + a.dual_d_y[aidx];
                slope[effective_k] = a.neg_dual_d_slope[aidx];
                effective_k++;
                aidx++;
            }
            aidx--; sumidx--;
            if (break_flag) break;
        } else { // a.neg_dual_d_slope[aidx] < sum_slope[sumidx]
            while (aidx < a.k && a.neg_dual_d_slope[aidx + 1] < sum_slope[sumidx]) {
                printf("%.10lf\n",a.neg_dual_d_slope[aidx + 1] - sum_slope[sumidx]);
                aidx++;
            }
            while (sumidx>=0 && a.neg_dual_d_slope[aidx + 1] >= sum_slope[sumidx]) {
                //a_i - a_pi = d_x[aidx], sum_x[sumidx] < a_i (< [sumidx+1])
                if (sum_x[sumidx] <= a.dual_d_x[aidx]){
                    dual_d_a_BT[effective_k] = (a.dual_d_x[aidx]);
                    slope[effective_k] = sum_slope[sumidx];
                    break_flag = true;
                    break;
                }
                dual_d_a_BT[effective_k] = sum_x[sumidx];
                x[effective_k] = sum_x[sumidx] - a.dual_d_x[aidx];
                y[effective_k] = sum_y[sumidx] + a.dual_d_y[aidx];
                slope[effective_k] = sum_slope[sumidx];
                effective_k++;
                sumidx--;
            }
            sumidx++; aidx++;
            if (break_flag) break;
        }
    }
    x[effective_k] = 0;
    y[effective_k] = y[effective_k-1] - x[effective_k-1] * slope[effective_k];
    std::reverse(x.begin(), x.begin() + effective_k+1);
    std::reverse(y.begin(), y.begin() + effective_k+1);
    std::reverse(slope.begin(), slope.begin() + effective_k+1);
    std::reverse(dual_d_a_BT.begin(), dual_d_a_BT.begin()+effective_k+1);
}

void DPSTATE::zero(PWL &a) {
    effective_k = std::upper_bound(a.dual_d_x.begin(),a.dual_d_x.end(),0) - a.dual_d_x.begin();
    if (effective_k>=a.k) y[0] = a.dual_d_y[effective_k-1] + a.dual_d_x[effective_k-1] * a.neg_dual_d_slope[effective_k];
    else y[0] = a.dual_d_y[effective_k] + a.dual_d_x[effective_k] * a.neg_dual_d_slope[effective_k];
    x[0] = 0;
    for (int i = 0; i < effective_k; i++){
        x[i+1] = -a.dual_d_x[effective_k - i - 1];
        y[i+1] = a.dual_d_y[effective_k - i - 1];
    }
    for (int i = 0; i <= effective_k; i++){
        slope[i] = a.neg_dual_d_slope[effective_k - i];
    }
}

real DPSTATE::backtrace(real val) {
    int idx = std::lower_bound(x.begin(), x.begin() + effective_k + 1, val) - x.begin();
    if (idx <= 0) return dual_d_a_BT[0];
    return (dual_d_a_BT[idx] - dual_d_a_BT[idx-1])/(x[idx]-x[idx-1])*(val-x[idx-1])+dual_d_a_BT[idx-1];
}

#undef dual_d_x
#undef neg_dual_d_slope