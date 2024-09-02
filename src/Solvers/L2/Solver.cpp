#include <algorithm>

#include "Solver.h"
#include "PiecewiseQuadraticF.h"
#include "../TreeStructuredDualDP.h"

namespace L2Solver {
    void Solver::solve() {
        size_t nrows = variant_reads.size();
        fs = forward_solve<PiecewiseQuadraticF>(clone_tree, vertex_map, variant_reads, total_reads, root);

        double obj = 0;
        for (size_t j = 0; j < nrows; ++j) {
            PiecewiseQuadraticF f = fs[vertex_map.at(root)][j];
            std::vector<double> cs = f.get_derivative_intercepts();

            double alpha_0 = 0.0;
            size_t i = 0;
            for (i = 0; i < cs.size();) {
                if (cs[i] - 1 <= 0) break; 
                i++;
            }

            if (i == 0) {
                // (c0 - 1) + m0 * alpha_0 = 0
                alpha_0 = (1.0 - f.c0) / f.m0;
            } else {
                // (c_i - 1) + m_i * (alpha_0 - b_i) = 0
                alpha_0 = (1.0 - cs[i - 1]) / f.slopes[i - 1] + f.breakpoints[i - 1];
            }

            alpha_0 = std::max(0.0, alpha_0);
            obj += f(alpha_0) - alpha_0;
            backtrack(alpha_0, j);
        }

        this->objective = obj;
    }

    void Solver::backtrack(double alpha_0, int j) {
        std::stack<int> stack;
        stack.push(root);

        while(!stack.empty()) {
            int i = stack.top();
            stack.pop();

            auto children = clone_tree.successors(i);

            PiecewiseQuadraticF g;
            for (auto k : clone_tree.successors(i)) {
                PiecewiseQuadraticF f = fs[vertex_map.at(k)][j];
                g = g + f;
                stack.push(k);
            }

            double gamma; 
            if (i == root) {
                gamma = alpha_0;
            } else {
                int p = clone_tree.predecessors(i)[0];
                gamma = alphas[j][p];
            }

            double freq = variant_reads[j][i] / total_reads[j][i];
            double alpha;
            if (children.size() == 0) {
                alpha = std::max(0.0, gamma - 2.0 * freq);
            } else {
                alpha = g.compute_argmin(gamma, variant_reads[j][i], total_reads[j][i]);
            }

            alphas[j][i] = alpha;
            frequencies[j][i] = (alpha - gamma) / 2.0 + freq;
        }
    }
};
