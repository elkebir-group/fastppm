#ifndef PIECEWISE_QUADRATICF_HPP
#define PIECEWISE_QUADRATICF_HPP

#define ADD_PRIORITY_QUEUE 0

#include <queue>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace L2Solver {
class PiecewiseQuadraticF {
public:
    // b_0 = -\infty, b_{k+1} = \infty, so are not stored
    float f0, c0, m0; 
    std::vector<float> breakpoints; // b_1 <= ... <= b_k
    std::vector<float> slopes;      // m_1 <= ... <= m_k

    PiecewiseQuadraticF(size_t n) {
        f0 = 0.0;
        c0 = 0.0;
        m0 = 0.0;
        breakpoints.reserve(n);
        slopes.reserve(n);
    }

    PiecewiseQuadraticF() = delete;

    void reset_to_leaf(float frequency, float weight) {
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
        
        float current_slope = m0;
        std::vector<size_t> idx(children.size(), 0);
        int added = 0;
        while (true) {
            float min_bp = std::numeric_limits<float>::infinity();
            for (size_t i = 0; i < children.size(); i++) {
                int c = children[i];
                if (idx[i] < fs[c].breakpoints.size()) {
                    float candidate = fs[c].breakpoints[idx[i]];
                    if (candidate < min_bp) {
                        min_bp = candidate;
                    }
                }
            }

            if (min_bp == std::numeric_limits<float>::infinity()) {
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

    void get_derivative_intercepts(std::vector<float>& cs_buffer) const {
        // compute intercepts of the pieces of the derivative, using continuity
        cs_buffer.resize(breakpoints.size());
        cs_buffer[0] = c0 + m0 * (breakpoints[0]);
        for (size_t i = 1; i < breakpoints.size(); i++) {
            cs_buffer[i] = cs_buffer[i - 1] + slopes[i - 1] * (breakpoints[i] - breakpoints[i - 1]);
        }
    }

    float operator()(float x, const std::vector<float>& cs) const {
        // first we find which piece x and 0.0 are within
        size_t k = 0, l = 0;
        for (size_t i = 0; i < breakpoints.size(); i++) {
            if (breakpoints[i] >= x) {
                k = i;
                break;
            }
        }

        for (size_t i = 0; i < breakpoints.size(); i++) {
            if (breakpoints[i] >= 0.0f) {
                l = i;
                break;
            }
        }

        // set value at the breakpoint containing 0.0
        float value;
        if (l == 0) {
            value = f0 + c0 * breakpoints[l] + 0.5f * m0 * (breakpoints[l] * breakpoints[l]);
        } else {
            value = f0 + (cs[l - 1] - slopes[l - 1] * breakpoints[l - 1]) * (breakpoints[l]) + 0.5 * slopes[l - 1] * (breakpoints[l] * breakpoints[l]);
        }

        // compute value at the breakpoint containing x
        if (k > l) {
            for (size_t i = l + 1; i < k; i++) {
                value += (cs[i - 1] - slopes[i - 1] * breakpoints[i - 1]) * (breakpoints[i] - breakpoints[i - 1]);
                value += 0.5f * slopes[i - 1] * (breakpoints[i] * breakpoints[i] - breakpoints[i - 1] * breakpoints[i - 1]);
            }
        }

        if (k <= l) {
            for (size_t i = l; i >= k && i > 0; i--) {
                value -= (cs[i - 1] - slopes[i - 1] * breakpoints[i - 1]) * (breakpoints[i] - breakpoints[i - 1]);
                value -= 0.5f * slopes[i - 1] * (breakpoints[i] * breakpoints[i] - breakpoints[i - 1] * breakpoints[i - 1]);
            }
        }

        // compute value at x
        if (k == 0) {
            return value + (c0 * (x - breakpoints[0]) + 0.5f * m0 * (x * x - breakpoints[0] * breakpoints[0]));
        } else {
            return value + (cs[k - 1] - slopes[k - 1] * breakpoints[k - 1]) * (x - breakpoints[k - 1]) + 0.5f * slopes[k - 1] * (x * x - breakpoints[k - 1] * breakpoints[k - 1]);
        }
    }

    // when F = \sum{j \in \delta(i)}J_j, this updates F to be 
    // J_i(\gamma) = max_{x \geq 0}(h_i(x - \gamma) + F(x))
    // really, this is the meat of the algorithm
    void update_representation(float frequency, float weight, PiecewiseQuadraticF& result, std::vector<float>& cs) const {
        // compute intercepts of the pieces of the derivative, using continuity
        float half_weight_inv = 1.0f / (2.0f * weight);

        get_derivative_intercepts(cs);

        // find first breakpoint x such that \alpha_i^*(x) = 0
        size_t l = 0;
        for (size_t i = 0; i < breakpoints.size(); i++) {
            if (breakpoints[i] > 0.0f) {
                l = i;
                break;
            }
        }

        float x = 0.0f;
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
            float new_breakpoint;
            if (i == 0) { 
                new_breakpoint = x;
            } else {
                new_breakpoint = 2.0f * weight * (frequency - cs[i + l - 1]) + breakpoints[i + l - 1];
            }

            float new_slope;
            if (i + l == 0) {
                new_slope = -m0 / (2.0f*weight*m0 - 1.0f);
            } else {
                new_slope = -(slopes[i + l - 1] / (2.0f*weight*slopes[i + l - 1] - 1.0f));
            }

            // remove duplicate slopes or breakpoints
            if (i > 0) {
                if ((new_breakpoint - result.breakpoints[i - num_duplicates - 1] < 1e-8f) && (new_breakpoint - result.breakpoints[i - num_duplicates - 1] > -1e-8f)) {
                    num_duplicates++;
                } else if ((new_slope - result.slopes[i - num_duplicates - 1] < 1e-9f) && (new_slope - result.slopes[i - num_duplicates - 1] > -1e-9f)) {
                    num_duplicates++;
                }
            }

            result.slopes[i - num_duplicates] = new_slope;
            result.breakpoints[i - num_duplicates] = new_breakpoint;
        }

        result.breakpoints.resize(result.breakpoints.size() - num_duplicates);
        result.slopes.resize(result.slopes.size() - num_duplicates);

        float new_m0 = -half_weight_inv;
        float new_c0 = frequency;

        for (size_t i = 0; i < breakpoints.size(); i++) {
            if (2.0 * weight * (frequency - cs[i]) + breakpoints[i] > 0.0f) {
                l = i;
                break;
            }
        }

        float alpha_star = 0.0f;

        if (l == 0) {
            alpha_star = (frequency - c0) / (m0 - half_weight_inv);
        } else {
            alpha_star = (frequency + half_weight_inv * breakpoints[l-1] - cs[l-1]) / (slopes[l-1] - half_weight_inv) + breakpoints[l-1];
        }

        alpha_star = std::max(0.0f, alpha_star);
        float new_f0 = this->operator()(alpha_star, cs) - (0.5f * half_weight_inv * alpha_star * alpha_star + frequency * alpha_star);

        result.f0 = new_f0;
        result.c0 = new_c0;
        result.m0 = new_m0;
    }

    float compute_argmin(float gamma, float frequency, float weight, std::vector<float> &cs) const {
        // compute intercepts of the pieces of the derivative, using continuity
        get_derivative_intercepts(cs);
        float half_weight_inv = 1.0 / (2.0 * weight);

        // find first breakpoint x
        size_t l = 0;
        for (;l < breakpoints.size(); l++) {
            if (2.0 * weight * (frequency - cs[l]) + breakpoints[l] > 0.0) {
                break;
            }
        }

        float alpha_star = 0.0;

        if (l == 0) {
            alpha_star = (frequency - half_weight_inv*gamma - c0) / (m0 - half_weight_inv);
        } else {
            alpha_star = (frequency - half_weight_inv*gamma + half_weight_inv * breakpoints[l-1] - cs[l-1]) / (slopes[l-1] - half_weight_inv) + breakpoints[l-1];
        }

        return std::max(0.0f, alpha_star);
    }
};
};

#endif
