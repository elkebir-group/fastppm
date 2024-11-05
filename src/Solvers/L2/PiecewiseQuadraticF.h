#ifndef PIECEWISE_QUADRATICF_HPP
#define PIECEWISE_QUADRATICF_HPP

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
    std::vector<float> slopes;      // m_1 <= ... <= m_k
    std::vector<float> breakpoints; // b_1 <= ... <= b_k

    PiecewiseQuadraticF() {
        f0 = 0.0;
        c0 = 0.0;
        m0 = 0.0;
    }

    // leaf constructor
    PiecewiseQuadraticF(float frequency, float weight) {
        breakpoints.push_back(2.0 * weight * frequency);
        slopes.push_back(0.0);
        f0 = 0.0;
        c0 = frequency;
        m0 = -1.0/(2.0*weight);
    }

    // runs in O(k + k') time
    PiecewiseQuadraticF operator+(const PiecewiseQuadraticF& other) const {
        PiecewiseQuadraticF result;
        result.f0 = f0 + other.f0;
        result.c0 = c0 + other.c0;
        result.m0 = m0 + other.m0;

        std::vector<float> merged_breakpoints(breakpoints.size() + other.breakpoints.size());
        std::vector<float> merged_slopes(slopes.size() + other.slopes.size());

        size_t l = 0, h = 0;
        for (size_t i = 0, j = 0; i < breakpoints.size() || j < other.breakpoints.size(); l++) {
            if (j == other.breakpoints.size() || (i != breakpoints.size() && breakpoints[i] < other.breakpoints[j])) {
                merged_breakpoints[l] = breakpoints[i];
                i++;
            } else if (i == breakpoints.size() || (j != other.breakpoints.size() && breakpoints[i] > other.breakpoints[j])) {
                merged_breakpoints[l] = other.breakpoints[j];
                j++;
            } else {
                if (i == breakpoints.size()) {
                    merged_breakpoints[l] = other.breakpoints[j];
                } else {
                    merged_breakpoints[l] = breakpoints[i];
                }

                i++;
                j++;
                h++;
            }

            float slope = 0;
            if (i == 0) {
                slope += m0;
            } else {
                slope += slopes[i - 1];
            }

            if (j == 0) {
                slope += other.m0;
            } else {
                slope += other.slopes[j - 1];
            }

            merged_slopes[l] = slope;
        }

        merged_breakpoints.resize(merged_breakpoints.size() - h);
        merged_slopes.resize(merged_slopes.size() - h);
        result.breakpoints = std::move(merged_breakpoints);
        result.slopes = std::move(merged_slopes);
        return result;
    }

    std::vector<float> get_derivative_intercepts() const {
        // compute intercepts of the pieces of the derivative, using continuity
        std::vector<float> cs(breakpoints.size());
        cs[0] = c0 + m0 * (breakpoints[0]);
        for (size_t i = 1; i < breakpoints.size(); i++) {
            cs[i] = cs[i - 1] + slopes[i - 1] * (breakpoints[i] - breakpoints[i - 1]);
        }

        return cs; 
    }

    float operator()(float x) const {
        // the above can all be done in O(k) time and precomputed
        std::vector<float> cs = get_derivative_intercepts();

        // first we find which piece x and 0.0 are within
        size_t k = std::lower_bound(breakpoints.begin(), breakpoints.end(), x) - breakpoints.begin();
        size_t l = std::lower_bound(breakpoints.begin(), breakpoints.end(), 0.0) - breakpoints.begin();

        // set value at the breakpoint containing 0.0
        float value;
        if (l == 0) {
            value = f0 + c0 * breakpoints[l] + 0.5 * m0 * (breakpoints[l] * breakpoints[l]);
        } else {
            value = f0 + (cs[l - 1] - slopes[l - 1] * breakpoints[l - 1]) * (breakpoints[l]) + 0.5 * slopes[l - 1] * (breakpoints[l] * breakpoints[l]);
        }

        // compute value at the breakpoint containing x
        if (k > l) {
            for (size_t i = l + 1; i < k; i++) {
                value += (cs[i - 1] - slopes[i - 1] * breakpoints[i - 1]) * (breakpoints[i] - breakpoints[i - 1]);
                value += 0.5 * slopes[i - 1] * (breakpoints[i] * breakpoints[i] - breakpoints[i - 1] * breakpoints[i - 1]);
            }
        }

        if (k <= l) {
            for (size_t i = l; i >= k && i > 0; i--) {
                value -= (cs[i - 1] - slopes[i - 1] * breakpoints[i - 1]) * (breakpoints[i] - breakpoints[i - 1]);
                value -= 0.5 * slopes[i - 1] * (breakpoints[i] * breakpoints[i] - breakpoints[i - 1] * breakpoints[i - 1]);
            }
        }

        // compute value at x
        if (k == 0) {
            return value + (c0 * (x - breakpoints[0]) + 0.5 * m0 * (x * x - breakpoints[0] * breakpoints[0]));
        } else {
            return value + (cs[k - 1] - slopes[k - 1] * breakpoints[k - 1]) * (x - breakpoints[k - 1]) + 0.5 * slopes[k - 1] * (x * x - breakpoints[k - 1] * breakpoints[k - 1]);
        }
    }

    // when F = \sum{j \in \delta(i)}J_j, this updates F to be 
    // J_i(\gamma) = max_{x \geq 0}(h_i(x - \gamma) + F(x))
    // really, this is the meat of the algorithm
    PiecewiseQuadraticF update_representation(float frequency, float weight) const {
        // compute intercepts of the pieces of the derivative, using continuity
        float half_weight_inv = 1.0 / (2.0 * weight);

        std::vector<float> cs = get_derivative_intercepts();

        // find first breakpoint x such that \alpha_i^*(x) = 0
        size_t l = std::lower_bound(breakpoints.begin(), breakpoints.end(), 0.0) - breakpoints.begin();
        float x = 0.0;
        if (l == 0) {
            x = frequency - c0;
        } else {
            x = frequency + (half_weight_inv * breakpoints[l-1]) - cs[l - 1] + (breakpoints[l-1] * (slopes[l-1] - half_weight_inv));
        }

        x *= 2.0 * weight;

        std::vector<float> zs(breakpoints);
        for (size_t i = 0; i < breakpoints.size(); i++) {
            zs[i] += 2.0 * weight * (frequency - cs[i]);
        }

        std::vector<float> new_breakpoints(breakpoints.size() + 1 - l);
        new_breakpoints[0] = x;
        for (size_t i = l; i < breakpoints.size(); i++) {
            new_breakpoints[i - l + 1] = zs[i];
        }

        std::vector<float> new_slopes(new_breakpoints.size());
        for (size_t i = 0; i < new_breakpoints.size(); i++) {
            float slope;
            if (i + l == 0) {
                slope = m0;
            } else {
                slope = slopes[i + l - 1];
            }

            new_slopes[i] = -(slope / (2*weight*slope - 1));
        }
        
        float new_m0 = -half_weight_inv;
        float new_c0 = frequency;

        l = std::lower_bound(zs.begin(), zs.end(), 0.0) - zs.begin();
        float alpha_star = 0.0;

        if (l == 0) {
            alpha_star = (frequency - c0) / (m0 - half_weight_inv);
        } else {
            alpha_star = (frequency + half_weight_inv * breakpoints[l-1] - cs[l-1]) / (slopes[l-1] - half_weight_inv) + breakpoints[l-1];
        }

        alpha_star = std::max(0.0f, alpha_star);
        float new_f0 = this->operator()(alpha_star) - (0.5 * half_weight_inv * alpha_star * alpha_star + frequency * alpha_star);

        std::vector<float> deduped_breakpoints(new_breakpoints.size());
        std::vector<float> deduped_slopes(new_slopes.size());
        size_t j = 0;
        for (size_t i = 0; i < new_breakpoints.size(); i++) {
            if (i == 0 || abs(new_breakpoints[i] - new_breakpoints[i - 1]) > 1e-8) { // to handle floating point errors
                deduped_breakpoints[j] = new_breakpoints[i];
                deduped_slopes[j] = new_slopes[i];
                j++;
            } else if (j > 0) {
                deduped_slopes[j - 1] = new_slopes[i];
            }
        }

        deduped_breakpoints.resize(j);
        deduped_slopes.resize(j);

        PiecewiseQuadraticF result;
        result.f0 = new_f0;
        result.c0 = new_c0;
        result.m0 = new_m0;
        result.breakpoints = deduped_breakpoints;
        result.slopes = deduped_slopes;

        return result;
    }

    float compute_argmin(float gamma, float frequency, float weight) const {
        // compute intercepts of the pieces of the derivative, using continuity
        std::vector<float> cs = get_derivative_intercepts();
        float half_weight_inv = 1.0 / (2.0 * weight);

        // find first breakpoint x
        std::vector<float> zs(breakpoints.size());
        for (size_t i = 0; i < breakpoints.size(); i++) {
            zs[i] = 2.0 * weight * (frequency - cs[i]) + breakpoints[i];
        }

        size_t l = std::lower_bound(zs.begin(), zs.end(), gamma) - zs.begin();
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
