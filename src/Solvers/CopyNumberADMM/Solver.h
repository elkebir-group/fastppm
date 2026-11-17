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
        // w_s, the weight of the squared loss of segment s on the state fractions; if empty, w_s = 1
        std::vector<double> segment_weight;
    };

    struct Options {
        double rho = 100.0;      // penalty, the same for every row of E
        // penalty of row j is kappa times the curvature of its loss; kappa < 0: kappa = 1 / sqrt(lambda_min lambda_max)
        // of Sigma0^(-1/2) E^T H E Sigma0^(-1/2) (root excluded); kappa = 0: the scalar penalty rho on every row
        double kappa = -1.0;
        double eps = 1.0;        // stop when L(u) - d(beta) <= eps (absolute, in units of the loss)
        int max_iterations = 100000;
        int check_every = 10;    // iterations between evaluations of the stopping rule
        int polish_steps = 1;    // exponentiated-gradient steps on the usages at each check (0: none)
        double alpha = 1.6;      // over-relaxation (Boyd et al. 2011, Sec. 3.4.3), in (0, 2); 1 is plain ADMM
        bool verbose = false;
        int anderson_memory = 10; // Anderson acceleration with this many past steps (0: off; SCS default 10)
        double anderson_reg = 1e-10;
        bool warm_prox = true;   // start the y-update of each SNV from its previous value
    };

    struct Result {
        std::vector<double> usages;       // u^(t), the usages of the dilution (25)
        std::vector<double> frequencies;  // their frequencies
        double loss;                      // L(u^(t))
        double lower_bound;               // d(beta^(t)) (24)
        int iterations;
        int anderson_accepted, anderson_rejected;
        double kappa;                     // kappa used (0 in scalar mode)
        bool converged;
        double runtime_ms;
    };

    Result solve(const Instance& instance, const Options& options);
};

#endif
