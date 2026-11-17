#ifndef COPY_NUMBER_ADMM_SOLVER_H
#define COPY_NUMBER_ADMM_SOLVER_H

#include <tuple>
#include <utility>
#include <vector>

/*
 * fastppm-CN: linearized ADMM (Algorithm 1 of the fastppm-CN paper) for
 * copy-number-aware clonal regression with the binomial loss on SNVs and the
 * squared loss on copy-number state fractions.
 *
 * Clones are 0..k-1. Equation numbers refer to the paper.
 */
namespace CopyNumberADMM {
    using State = std::pair<int, int>;  // allele-specific copy numbers (a, b)

    struct Instance {
        std::vector<int> parent;                       // parent[v], -1 at the root
        std::vector<std::vector<State>> states;        // states[s][v] = (a(v,s), b(v,s))
        std::vector<std::vector<int>> multiplicity;    // multiplicity[i][v] = m(v,i)
        std::vector<int> segment;                      // seg(i)
        std::vector<int> variant, total;               // v_i, d_i
        // copy-number calls: per segment, (a, b, mu_{s,(a,b)}); states of C_s not listed get 0
        std::vector<std::vector<std::tuple<int, int, double>>> state_fractions;
        // hat c_s; if empty, hat c_s = sum over (a,b) of (a + b) mu_{s,(a,b)}
        std::vector<double> mean_copy_number;
    };

    struct Options {
        double rho = 100.0;      // penalty, the same for every row of E
        double kappa = 0.0;      // if > 0, penalty of row j is kappa times the curvature of its loss
        double eps = 1.0;        // stop when L(u) - d(beta) <= eps (absolute, in units of the loss)
        int max_iterations = 100000;
        int check_every = 10;    // iterations between evaluations of the stopping rule
        bool verbose = false;
    };

    struct Result {
        std::vector<double> usages;       // u^(t), the usages of the dilution (25)
        std::vector<double> frequencies;  // their frequencies
        double loss;                      // L(u^(t))
        double lower_bound;               // d(beta^(t)) (24)
        int iterations;
        bool converged;
        double runtime_ms;
    };

    Result solve(const Instance& instance, const Options& options);
};

#endif
