#ifndef PIECEWISE_QUADRATICF_HPP
#define PIECEWISE_QUADRATICF_HPP

#define ADD_PRIORITY_QUEUE 0

#include <queue>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

namespace L2Solver {
template <typename T>
class PiecewiseQuadraticF {
public:
    // b_0 = -\infty, b_{k+1} = \infty, so are not stored
    T f0, c0, m0; 
    std::vector<T> breakpoints; // b_1 <= ... <= b_k
    std::vector<T> slopes;      // m_1 <= ... <= m_k

    PiecewiseQuadraticF(size_t n) {
        f0 = 0.0;
        c0 = 0.0;
        m0 = 0.0;
        breakpoints.reserve(n);
        slopes.reserve(n);
    }

    PiecewiseQuadraticF() = delete;

    void reset_to_leaf(T frequency, T weight) {
        f0 = 0.0;
        c0 = frequency;
        m0 = -1.0 / (2.0 * weight);
        breakpoints.clear();
        breakpoints.push_back(2.0f * weight * frequency);
        slopes.clear();
        slopes.push_back(0.0f);
    }

    void sum(const std::vector<int>& children, const std::vector<PiecewiseQuadraticF>& fs) {
        f0 = 0.0f;
        c0 = 0.0f;
        m0 = 0.0f;
        breakpoints.clear();
        slopes.clear();

        if (children.empty()) return;

        if (children.size() == 1) {
            f0 = fs[children[0]].f0;
            c0 = fs[children[0]].c0;
            m0 = fs[children[0]].m0;
            breakpoints = fs[children[0]].breakpoints;
            slopes = fs[children[0]].slopes;
            return;
        }
        
        for (int c : children) {
            f0 += fs[c].f0;
            c0 += fs[c].c0;
            m0 += fs[c].m0;
        }

        size_t total_bp_count = 0;
        for (int c : children) {
          total_bp_count += fs[c].breakpoints.size();
        }

        breakpoints.resize(total_bp_count);
        slopes.resize(total_bp_count);
        
        T current_slope = m0;
        std::vector<size_t> idx(children.size(), 0);
        int added = 0;
        while (true) {
            T min_bp = std::numeric_limits<T>::infinity();
            for (size_t i = 0; i < children.size(); i++) {
                int c = children[i];
                if (idx[i] < fs[c].breakpoints.size()) {
                    T candidate = fs[c].breakpoints[idx[i]];
                    if (candidate < min_bp) {
                        min_bp = candidate;
                    }
                }
            }

            if (min_bp == std::numeric_limits<T>::infinity()) {
                break;
            }

            for (size_t i = 0; i < children.size(); i++) {
                int c = children[i];
                if (idx[i] < fs[c].breakpoints.size() && fs[c].breakpoints[idx[i]] == min_bp) {
                    if (idx[i] == 0) {
                        current_slope -= fs[c].m0;
                    } else {
                        current_slope -= fs[c].slopes[idx[i] - 1];
                    }
                    current_slope += fs[c].slopes[idx[i]];
                    idx[i]++;
                }
            }
            
            breakpoints[added] = min_bp;
            slopes[added] = current_slope;
            added++;
        }

        breakpoints.resize(added);
        slopes.resize(added);
    }

    void get_derivative_intercepts(std::vector<T>& cs_buffer) const {
        // compute intercepts of the pieces of the derivative, using continuity
        cs_buffer.resize(breakpoints.size());
        cs_buffer[0] = c0 + m0 * (breakpoints[0]);
        for (size_t i = 1; i < breakpoints.size(); i++) {
            cs_buffer[i] = cs_buffer[i - 1] + slopes[i - 1] * (breakpoints[i] - breakpoints[i - 1]);
        }
    }

    // the derivative on piece p, which is (-inf, b_1) for p = 0 and [b_p, b_{p+1}) otherwise
    T derivative(T x, size_t p, const std::vector<T>& cs) const {
        return p == 0 ? c0 + m0 * x : cs[p - 1] + slopes[p - 1] * (x - breakpoints[p - 1]);
    }

    T slope(size_t p) const { return p == 0 ? m0 : slopes[p - 1]; }

    // F(x) = F(0) + the integral of F' from 0 to x, where F(0) = f0
    T operator()(T x, const std::vector<T>& cs) const {
        T lo = std::min(T(0), x), hi = std::max(T(0), x);
        size_t p = std::upper_bound(breakpoints.begin(), breakpoints.end(), lo) - breakpoints.begin();
        T integral = 0, t = lo;
        while (t < hi) {
            T next = p < breakpoints.size() ? std::min(breakpoints[p], hi) : hi;
            integral += derivative(t, p, cs) * (next - t) + T(0.5) * slope(p) * (next - t) * (next - t);
            t = next;
            p++;
        }
        return x >= 0 ? f0 + integral : f0 - integral;
    }

    // when F = \sum{j \in \delta(i)}J_j, this updates F to be 
    // J_i(\gamma) = max_{x \geq 0}(h_i(x - \gamma) + F(x))
    // really, this is the meat of the algorithm
    void update_representation(T frequency, T weight, PiecewiseQuadraticF& result, std::vector<T>& cs) const {
        // compute intercepts of the pieces of the derivative, using continuity
        T half_weight_inv = 1.0f / (2.0f * weight);

        get_derivative_intercepts(cs);

        // find the piece of F that contains 0: l is the number of breakpoints <= 0
        size_t l = breakpoints.size();
        for (size_t i = 0; i < breakpoints.size(); i++) {
            if (breakpoints[i] > 0.0f) {
                l = i;
                break;
            }
        }

        T x = 0.0f;
        if (l == 0) {
            x = frequency - c0;
        } else {
            x = frequency + (half_weight_inv * breakpoints[l-1]) - cs[l - 1] + (breakpoints[l-1] * (slopes[l-1] - half_weight_inv));
        }

        x *= 2.0f * weight;

        result.breakpoints.resize(breakpoints.size() + 1 - l);
        result.slopes.resize(result.breakpoints.size());
        int num_duplicates = 0;
        for (size_t i = 0; i < result.breakpoints.size(); i++) {
            T new_breakpoint;
            if (i == 0) { 
                new_breakpoint = x;
            } else {
                new_breakpoint = 2.0f * weight * (frequency - cs[i + l - 1]) + breakpoints[i + l - 1];
            }

            T new_slope;
            if (i + l == 0) {
                new_slope = -m0 / (2.0f*weight*m0 - 1.0f);
            } else {
                new_slope = -(slopes[i + l - 1] / (2.0f*weight*slopes[i + l - 1] - 1.0f));
            }

            // remove duplicate breakpoints (the earlier piece has no width, so the new piece replaces it)
            // and duplicate slopes (the new breakpoint is not a kink, so it is dropped)
            // (tolerances relative to the magnitudes, which vary with the weights)
            if (i > 0) {
                T prev_breakpoint = result.breakpoints[i - num_duplicates - 1];
                T prev_slope = result.slopes[i - num_duplicates - 1];
                T tol = 64 * std::numeric_limits<T>::epsilon();
                if (std::abs(new_breakpoint - prev_breakpoint) <= tol * std::max(std::abs(new_breakpoint), std::abs(prev_breakpoint))) {
                    num_duplicates++;
                } else if (std::abs(new_slope - prev_slope) <= tol * std::max(std::abs(new_slope), std::abs(prev_slope))) {
                    num_duplicates++;
                    continue;
                }
            }

            result.slopes[i - num_duplicates] = new_slope;
            result.breakpoints[i - num_duplicates] = new_breakpoint;
        }

        result.breakpoints.resize(result.breakpoints.size() - num_duplicates);
        result.slopes.resize(result.slopes.size() - num_duplicates);

        T new_m0 = -half_weight_inv;
        T new_c0 = frequency;

        // find the piece of F that contains alpha^*(0): l is the number of breakpoints of the new function <= 0
        l = breakpoints.size();
        for (size_t i = 0; i < breakpoints.size(); i++) {
            if (2.0 * weight * (frequency - cs[i]) + breakpoints[i] > 0.0f) {
                l = i;
                break;
            }
        }

        T alpha_star = 0.0f;

        if (l == 0) {
            alpha_star = (frequency - c0) / (m0 - half_weight_inv);
        } else {
            alpha_star = (frequency + half_weight_inv * breakpoints[l-1] - cs[l-1]) / (slopes[l-1] - half_weight_inv) + breakpoints[l-1];
        }

        alpha_star = std::max(T(0), alpha_star);
        T new_f0 = this->operator()(alpha_star, cs) - (0.5f * half_weight_inv * alpha_star * alpha_star + frequency * alpha_star);

        result.f0 = new_f0;
        result.c0 = new_c0;
        result.m0 = new_m0;
    }

    T compute_argmin(T gamma, T frequency, T weight, std::vector<T> &cs) const {
        // compute intercepts of the pieces of the derivative, using continuity
        get_derivative_intercepts(cs);
        T half_weight_inv = 1.0 / (2.0 * weight);

        // find first breakpoint x
        size_t l = 0;
        for (;l < breakpoints.size(); l++) {
            if (2.0 * weight * (frequency - cs[l]) + breakpoints[l] >= gamma) {
                break;
            }
        }

        T alpha_star = 0.0;

        if (l == 0) {
            alpha_star = (frequency - half_weight_inv*gamma - c0) / (m0 - half_weight_inv);
        } else {
            alpha_star = (frequency - half_weight_inv*gamma + half_weight_inv * breakpoints[l-1] - cs[l-1]) / (slopes[l-1] - half_weight_inv) + breakpoints[l-1];
        }

        return std::max(T(0), alpha_star);
    }
};
};

#endif
