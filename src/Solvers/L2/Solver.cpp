#include <algorithm>

#include "Solver.h"
#include "PiecewiseQuadraticF.h"
#include "../TreeStructuredDualDP.h"

namespace L2Solver {
    template <typename T>
    void BasicSolver<T>::initialize() {
        postorder = clone_tree.postorder_traversal(vertex_map[root]);
        preorder = clone_tree.preorder_traversal(vertex_map[root]);
        num_descendants = std::vector(clone_tree.size(), 1);
        for (auto u : postorder) {
            if (clone_tree.out_degree(u) == 0) continue;

            const auto& children = clone_tree.successors(u);
            for (auto v : children) {
                num_descendants[u] += num_descendants[v];
            }
        }

        for (size_t i = 0; i < clone_tree.size(); ++i) {
            fs.push_back(PiecewiseQuadraticF<T>(num_descendants[i] + 1));
            gs.push_back(PiecewiseQuadraticF<T>(num_descendants[i] + 1));
        }

        cs.reserve(clone_tree.size() + 1);
    }

    template <typename T>
    void BasicSolver<T>::solve() {
        size_t nrows = frequency_matrix.size();

        T obj = 0;
        for (size_t j = 0; j < nrows; ++j) {
            forward_solve<PiecewiseQuadraticF<T>, T>(clone_tree, postorder, vertex_map, frequency_matrix, weight_matrix, root, fs, gs, cs, j);
            const PiecewiseQuadraticF<T>& f = fs[vertex_map[root]];
            f.get_derivative_intercepts(cs);

            T alpha_0 = 0.0;
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

            // f_root <= 1 needs alpha_0 >= 0; f_root = 1 leaves alpha_0 free
            if (!fix_root_frequency) alpha_0 = std::max(T(0), alpha_0);
            obj += f(alpha_0, cs) - alpha_0;
            backtrack(alpha_0, j);
        }

        this->objective = obj;
    }

    template <typename T>
    void BasicSolver<T>::backtrack(T alpha_0, int j) {
      for (auto u : preorder) {
            int i = clone_tree[u].data;

            const auto& children = clone_tree.successors(u);

            const PiecewiseQuadraticF<T>& g = gs[u];

            T gamma;
            if (i == root) {
                gamma = alpha_0;
            } else {
                int p = clone_tree.predecessors(u)[0]; // convert i to vertex coordinates, then get the parent
                gamma = alphas[j][p]; // get the parent's alpha
            }

            T freq = frequency_matrix[j][i];
            T weight = weight_matrix[j][i];

            T alpha;
            if (children.size() == 0) {
                alpha = std::max(T(0), gamma - T(2) * weight * freq);
            } else {
                alpha = g.compute_argmin(gamma, freq, weight, cs);
            }

            alphas[j][u] = alpha;
            frequencies[j][i] = (alpha - gamma) / (2.0 * weight) + freq;
        }
    }

    template class BasicSolver<float>;
    template class BasicSolver<double>;
};
