#include <algorithm>

#include "Solver.h"
#include "PiecewiseQuadraticF.h"
#include "../TreeStructuredDualDP.h"

namespace L2Solver {
    void Solver::solve() {
        size_t nrows = frequency_matrix.size();
        
        forward_solve<PiecewiseQuadraticF>(clone_tree, vertex_map, frequency_matrix, weight_matrix, root, fs, gs);

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
                alpha_0 = (1.0 - f.c0) / f.m0;
            } else {
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
        stack.push(root); // root is in column coordinates

        while(!stack.empty()) {
            int i = stack.top(); // i is in column coordinates
            stack.pop();

            auto children = clone_tree.successors(vertex_map.at(i));

            PiecewiseQuadraticF g = gs[i][j];

            double gamma; 
            if (i == root) {
                gamma = alpha_0;
            } else {
                int p = clone_tree.predecessors(vertex_map.at(i))[0]; // convert i to vertex coordinates, then get the parent
                gamma = alphas[j][clone_tree[p].data]; // get the parent's alpha
            }

            double freq = frequency_matrix[j][i];
            double weight = weight_matrix[j][i];

            double alpha;
            if (children.size() == 0) {
                alpha = std::max(0.0, gamma - 2.0 * weight * freq);
            } else {
                alpha = g.compute_argmin(gamma, freq, weight);
            }

            alphas[j][i] = alpha;
            frequencies[j][i] = (alpha - gamma) / (2.0 * weight) + freq;
        }
    }
};
