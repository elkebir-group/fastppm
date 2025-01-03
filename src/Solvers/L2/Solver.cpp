#include <algorithm>

#include "Solver.h"
#include "PiecewiseQuadraticF.h"
#include "../TreeStructuredDualDP.h"

namespace L2Solver {
    void Solver::solve() {
        size_t nrows = frequency_matrix.size();
        size_t ncols = frequency_matrix[0].size();

        for (size_t i = 0; i < ncols; ++i) {
            fs.push_back(PiecewiseQuadraticF(ncols + 1));
            gs.push_back(PiecewiseQuadraticF(ncols + 1));
        }

        std::vector<int> postorder = clone_tree.postorder_traversal(vertex_map[root]);
        std::vector<int> preorder = clone_tree.preorder_traversal(vertex_map[root]);
        
        float obj = 0;
        for (size_t j = 0; j < nrows; ++j) {
            forward_solve<PiecewiseQuadraticF>(clone_tree, postorder, vertex_map, frequency_matrix, weight_matrix, root, fs, gs, j);
            const PiecewiseQuadraticF& f = fs[vertex_map[root]];
            const std::vector<float> cs = f.get_derivative_intercepts();

            float alpha_0 = 0.0;
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

            alpha_0 = std::max(0.0f, alpha_0);
            obj += f(alpha_0, cs) - alpha_0;
            backtrack(preorder, alpha_0, j);
        }

        this->objective = obj;
    }

    void Solver::backtrack(const std::vector<int>& preorder, float alpha_0, int j) {
      for (auto u : preorder) {
            int i = clone_tree[u].data;

            const auto& children = clone_tree.successors(u);

            const PiecewiseQuadraticF& g = gs[u];

            float gamma; 
            if (i == root) {
                gamma = alpha_0;
            } else {
                int p = clone_tree.predecessors(u)[0]; // convert i to vertex coordinates, then get the parent
                gamma = alphas[j][clone_tree[p].data]; // get the parent's alpha
            }

            float freq = frequency_matrix[j][i];
            float weight = weight_matrix[j][i];

            float alpha;
            if (children.size() == 0) {
                alpha = std::max(0.0f, gamma - 2.0f * weight * freq);
            } else {
                alpha = g.compute_argmin(gamma, freq, weight);
            }

            alphas[j][i] = alpha;
            frequencies[j][i] = (alpha - gamma) / (2.0 * weight) + freq;
        }
    }
};
