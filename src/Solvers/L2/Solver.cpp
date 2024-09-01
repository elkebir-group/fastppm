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
            std::vector<double> breakpoints = f.breakpoints;
            std::vector<double> values = f.evaluate_at_breakpoints();

            double alpha_0 = 0.0;
            double obj_inc = f.f0;
            for (size_t i = 0; i < breakpoints.size(); ++i) {
                if (breakpoints[i] < 0.0) continue;
                if (obj_inc < values[i] - breakpoints[i]) {
                    alpha_0 = breakpoints[i];
                    obj_inc = values[i] - breakpoints[i];
                }
            }

            obj += obj_inc;
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

            double freq = (double)variant_reads[j][i] / (double)total_reads[j][i];
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
