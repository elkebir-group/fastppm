#include "Solver.h"
#include "PiecewiseQuadraticF.h"
#include "../TreeStructuredDualDP.h"

namespace L2Solver {
    void Solver::solve() {
        size_t nrows = variant_reads.size();
        fs = forward_solve<PiecewiseQuadraticF>(clone_tree, vertex_map, variant_reads, total_reads, root);

        float obj = 0;
        for (size_t j = 0; j < nrows; ++j) {
            PiecewiseQuadraticF f = fs[vertex_map.at(root)][j];
            std::vector<double> breakpoints = f.breakpoints;
            std::vector<double> values = f.evaluate_at_breakpoints();
            double obj_inc = f.f0;
            for (size_t i = 0; i < breakpoints.size(); ++i) {
                if (breakpoints[i] < 0.0) continue;
                obj_inc = std::max(obj_inc, values[i] - breakpoints[i]);
            }

            obj += obj_inc;
        }

        this->objective = obj;
    }
};
