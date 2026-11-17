#include "Solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>

#include <spdlog/spdlog.h>

#include "DiGraph.h"
#include "../L2/Solver.h"
#include "../../Polynomial.h"

namespace CopyNumberADMM {
namespace {
    const double INF = std::numeric_limits<double>::infinity();

    // Compressed sparse rows.
    struct Sparse {
        int rows = 0, cols = 0;
        std::vector<int> ptr, idx;
        std::vector<double> val;

        void multiply(const std::vector<double>& x, std::vector<double>& y) const {
            y.assign(rows, 0.0);
            for (int r = 0; r < rows; ++r) {
                double s = 0.0;
                for (int p = ptr[r]; p < ptr[r + 1]; ++p) s += val[p] * x[idx[p]];
                y[r] = s;
            }
        }

        Sparse transpose() const {
            Sparse t;
            t.rows = cols;
            t.cols = rows;
            t.ptr.assign(cols + 1, 0);
            for (int c : idx) t.ptr[c + 1]++;
            for (int c = 0; c < cols; ++c) t.ptr[c + 1] += t.ptr[c];
            t.idx.resize(idx.size());
            t.val.resize(val.size());
            std::vector<int> next(t.ptr.begin(), t.ptr.end() - 1);
            for (int r = 0; r < rows; ++r) {
                for (int p = ptr[r]; p < ptr[r + 1]; ++p) {
                    int q = next[idx[p]]++;
                    t.idx[q] = r;
                    t.val[q] = val[p];
                }
            }
            return t;
        }
    };

    double xlogy(double x, double y) { return x == 0.0 ? 0.0 : x * std::log(y); }

    // Binomial negative log-likelihood of the VAF p, +inf outside [0, 1].
    double binomial_loss(double p, int v, int d) {
        if (p < 0.0 || p > 1.0) return INF;
        return -xlogy(v, p) - xlogy(d - v, 1.0 - p);
    }

    // argmin over p in [0, 1] of binomial_loss(p) + r/2 (p - w)^2, i.e. the
    // y-update (12) in the scaled variable p = y / hat c with r = rho hat c^2.
    // Safeguarded Newton steps on the derivative, which is increasing on (0, 1),
    // started from `start` if it lies in (0, 1) (the previous iterate) and
    // otherwise from a root of the cubic (16) by Cardano's formula.
    double binomial_prox(double r, double w, int v, int d, double start) {
        if (d == 0) return std::clamp(w, 0.0, 1.0);
        if (v == 0 && d - r * w >= 0.0) return 0.0;
        if (v == d && r * (1.0 - w) - d <= 0.0) return 1.0;

        double p = 0.5;
        if (start > 0.0 && start < 1.0) {
            p = start;
        } else {
            double x1 = -1, x2 = -1, x3 = -1;
            if (v > 0 && v < d) {
                solve_cubic(-r, r * (1.0 + w), d - r * w, -static_cast<double>(v), x1, x2, x3);
            } else if (v == 0) {
                solve_quadratic(-r, r * (1.0 + w), d - r * w, x1, x2);
            } else {
                solve_quadratic(r, -r * w, -static_cast<double>(d), x1, x2);
            }
            for (double x : {x1, x2, x3}) {
                if (x > 0.0 && x < 1.0) { p = x; break; }
            }
        }

        double lo = 0.0, hi = 1.0;
        for (int it = 0; it < 100; ++it) {
            double phi = r * (p - w) - v / p + (d - v) / (1.0 - p);
            if (phi > 0.0) hi = p; else lo = p;
            double dphi = r + v / (p * p) + (d - v) / ((1.0 - p) * (1.0 - p));
            double next = p - phi / dphi;
            if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
            if (std::abs(next - p) <= 1e-15 * std::max(1.0, p)) { p = next; break; }
            p = next;
        }
        return p;
    }

    // h_i(beta) = min over y in [0, c] of L_i(y) - beta y.
    double h_snv(double beta, double c, int v, int d) {
        double y;
        if (d == 0) {
            y = beta > 0.0 ? c : 0.0;
        } else if (v == 0) {
            y = beta <= d / c ? 0.0 : c - d / beta;
        } else if (v == d) {
            y = beta >= -d / c ? c : -d / beta;
        } else if (beta == 0.0) {
            y = v * c / d;
        } else {
            // stationary point: beta y^2 + (d - beta c) y - v c = 0, a single root in (0, c)
            double a = beta, b = d - beta * c, cc = -v * c;
            double disc = std::sqrt(std::max(b * b - 4.0 * a * cc, 0.0));
            double q = -0.5 * (b + (b >= 0.0 ? disc : -disc));
            double r1 = q / a, r2 = q != 0.0 ? cc / q : 0.0;
            y = a > 0.0 ? std::max(r1, r2) : std::min(r1, r2);
            y = std::clamp(y, 0.0, c);
        }
        return binomial_loss(y / c, v, d) - beta * y;
    }

    // Euclidean projection onto the probability simplex (Duchi et al. 2008).
    void project_simplex(std::vector<double>& w) {
        std::vector<double> s(w);
        std::sort(s.begin(), s.end(), std::greater<double>());
        double css = 0.0, tau = 0.0;
        for (size_t j = 0; j < s.size(); ++j) {
            css += s[j];
            double t = (css - 1.0) / static_cast<double>(j + 1);
            if (s[j] - t > 0.0) tau = t;
        }
        for (double& x : w) x = std::max(x - tau, 0.0);
    }
}

Result solve(const Instance& in, const Options& opt) {
    auto start = std::chrono::high_resolution_clock::now();
    const int k = static_cast<int>(in.parent.size());
    const int n = static_cast<int>(in.multiplicity.size());
    const int ell = static_cast<int>(in.states.size());

    int root = -1;
    std::vector<std::vector<int>> children(k);
    for (int v = 0; v < k; ++v) {
        if (in.parent[v] < 0) root = v; else children[in.parent[v]].push_back(v);
    }
    if (root < 0) throw std::invalid_argument("no root");
    std::vector<int> preorder{root};
    for (size_t j = 0; j < preorder.size(); ++j) {
        for (int w : children[preorder[j]]) preorder.push_back(w);
    }
    if (static_cast<int>(preorder.size()) != k) throw std::invalid_argument("parent does not define a tree");

    // Event matrix E (Definition 2): SNV rows, then one row per (s, (a, b)).
    Sparse E;
    E.cols = k;
    E.ptr.push_back(0);
    auto push_row = [&](const std::vector<double>& row) {
        for (int v = 0; v < k; ++v) {
            if (row[v] != 0.0) { E.idx.push_back(v); E.val.push_back(row[v]); }
        }
        E.ptr.push_back(static_cast<int>(E.idx.size()));
        E.rows++;
    };
    std::vector<double> row(k);
    for (int i = 0; i < n; ++i) {
        for (int v = 0; v < k; ++v) {
            int par = in.parent[v] < 0 ? 0 : in.multiplicity[i][in.parent[v]];
            row[v] = in.multiplicity[i][v] - par;
        }
        push_row(row);
    }
    std::vector<int> seg_start(ell + 1, n);
    std::vector<double> mu;                       // mu_{s,(a,b)} in row order
    std::vector<double> chat_seg(ell, 0.0);
    for (int s = 0; s < ell; ++s) {
        std::map<State, double> calls;
        if (s < static_cast<int>(in.state_fractions.size())) {
            for (auto& [a, b, m] : in.state_fractions[s]) calls[{a, b}] += m;
        }
        std::map<State, int> Cs;                  // C_s, in sorted order
        for (const State& st : in.states[s]) Cs[st] = 0;
        for (auto& [st, unused] : Cs) {
            for (int v = 0; v < k; ++v) {
                double chi = in.states[s][v] == st ? 1.0 : 0.0;
                double par = in.parent[v] < 0 ? 0.0 : (in.states[s][in.parent[v]] == st ? 1.0 : 0.0);
                row[v] = chi - par;
            }
            push_row(row);
            double m = calls.count(st) ? calls[st] : 0.0;
            mu.push_back(m);
            chat_seg[s] += (st.first + st.second) * m;
        }
        seg_start[s + 1] = E.rows;
    }
    if (!in.mean_copy_number.empty()) chat_seg = in.mean_copy_number;
    std::vector<double> chat(n);
    for (int i = 0; i < n; ++i) {
        chat[i] = chat_seg[in.segment[i]];
        if (!(chat[i] > 0.0)) throw std::invalid_argument("mean copy number of segment " + std::to_string(in.segment[i]) + " is not positive");
    }
    const Sparse ET = E.transpose();
    const int rows = E.rows;

    // Penalty per row of E: rho, or kappa times the curvature of the row's loss
    // (SNV i: d / (hat c^2 p (1 - p)) at p = (v + 1/2) / (d + 1); segment rows: 2).
    // A diagonal penalty R is ADMM on the rescaled constraint R^(1/2)(y, g) = R^(1/2) E f.
    std::vector<double> rho(rows, opt.rho);
    if (opt.kappa > 0.0) {
        for (int i = 0; i < n; ++i) {
            double p = (in.variant[i] + 0.5) / (in.total[i] + 1.0);
            rho[i] = opt.kappa * std::max(in.total[i], 1) / (chat[i] * chat[i] * p * (1.0 - p));
        }
        for (int r = n; r < rows; ++r) rho[r] = opt.kappa * 2.0;
    }

    // Damping (19), sigma_v = sum_j rho_j |e_jv| sum_w |e_jw|.
    std::vector<double> row_abs(rows, 0.0), sigma(k, 0.0);
    for (int r = 0; r < rows; ++r) {
        for (int p = E.ptr[r]; p < E.ptr[r + 1]; ++p) row_abs[r] += std::abs(E.val[p]);
    }
    for (int v = 0; v < k; ++v) {
        for (int p = ET.ptr[v]; p < ET.ptr[v + 1]; ++p) sigma[v] += rho[ET.idx[p]] * std::abs(ET.val[p]) * row_abs[ET.idx[p]];
        if (!(sigma[v] > 0.0)) throw std::invalid_argument("clone " + std::to_string(v) + " does not differ from its parent");
    }

    // f-update (20) by TSDDP with weights sigma/2 and the root frequency fixed at 1.
    digraph<int> tree;
    std::vector<int> vertex_map(k);
    for (int v = 0; v < k; ++v) vertex_map[v] = tree.add_vertex(v);
    for (int v = 0; v < k; ++v) {
        if (in.parent[v] >= 0) tree.add_edge(vertex_map[in.parent[v]], vertex_map[v]);
    }
    std::vector<std::vector<double>> target(1, std::vector<double>(k, 0.0));
    std::vector<std::vector<double>> weight(1, std::vector<double>(k));
    for (int v = 0; v < k; ++v) weight[0][v] = sigma[v] / 2.0;
    L2Solver::BasicSolver<double> tsddp(tree, vertex_map, target, weight, root);
    tsddp.set_fix_root_frequency(true);
    tsddp.initialize();

    auto f_from_u = [&](const std::vector<double>& u) {
        std::vector<double> f(u);
        for (int j = k - 1; j > 0; --j) f[in.parent[preorder[j]]] += f[preorder[j]];
        return f;
    };
    auto u_from_f = [&](const std::vector<double>& f) {
        std::vector<double> u(f);
        for (int v = 0; v < k; ++v) {
            if (in.parent[v] >= 0) u[in.parent[v]] -= f[v];
        }
        return u;
    };

    // Total loss (3) of the usages of f.
    std::vector<double> Ef_buf(rows);
    auto loss = [&](const std::vector<double>& f) {
        E.multiply(f, Ef_buf);
        double L = 0.0;
        for (int i = 0; i < n; ++i) L += binomial_loss(Ef_buf[i] / chat[i], in.variant[i], in.total[i]);
        for (int r = n; r < rows; ++r) L += (Ef_buf[r] - mu[r - n]) * (Ef_buf[r] - mu[r - n]);
        return L;
    };

    // Lower bound d(beta) (24).
    std::vector<double> psi(k), path(k), gs;
    auto dual_value = [&](const std::vector<double>& beta) {
        double d = 0.0;
        for (int i = 0; i < n; ++i) d += h_snv(beta[i], chat[i], in.variant[i], in.total[i]);
        for (int s = 0; s < ell; ++s) {
            gs.assign(seg_start[s + 1] - seg_start[s], 0.0);
            for (int r = seg_start[s]; r < seg_start[s + 1]; ++r) gs[r - seg_start[s]] = mu[r - n] + beta[r] / 2.0;
            project_simplex(gs);
            for (int r = seg_start[s]; r < seg_start[s + 1]; ++r) {
                double g = gs[r - seg_start[s]];
                d += (g - mu[r - n]) * (g - mu[r - n]) - beta[r] * g;
            }
        }
        ET.multiply(beta, psi);
        double best = INF;
        for (int v : preorder) {
            path[v] = psi[v] + (in.parent[v] < 0 ? 0.0 : path[in.parent[v]]);
            best = std::min(best, path[v]);
        }
        return d + best;
    };

    // Dilution (25).
    auto dilute = [&](const std::vector<double>& f) {
        E.multiply(f, Ef_buf);
        double theta = 1.0;
        for (int i = 0; i < n; ++i) {
            if (Ef_buf[i] > 0.0) theta = std::min(theta, chat[i] / Ef_buf[i]);
        }
        std::vector<double> fd(f);
        for (double& x : fd) x *= theta;
        fd[root] += 1.0 - theta;
        return fd;
    };

    // Algorithm 1.
    std::vector<double> f = f_from_u(std::vector<double>(k, 1.0 / k));
    std::vector<double> lam(rows, 0.0), Ef(rows), z(rows), yg(rows), q(rows), ETq(k), beta(rows), w;
    E.multiply(f, Ef);

    Result res;
    res.converged = false;
    std::vector<double> fd;
    double L = INF, d = -INF;
    int t = 0;
    while (t < opt.max_iterations) {
        for (int r = 0; r < rows; ++r) z[r] = Ef[r] + lam[r];
        for (int i = 0; i < n; ++i) {
            yg[i] = chat[i] * binomial_prox(rho[i] * chat[i] * chat[i], z[i] / chat[i], in.variant[i], in.total[i], -1.0);
        }
        for (int s = 0; s < ell; ++s) {
            w.assign(seg_start[s + 1] - seg_start[s], 0.0);
            const double rs = rho[seg_start[s]];  // constant within a segment
            for (int r = seg_start[s]; r < seg_start[s + 1]; ++r) w[r - seg_start[s]] = (2.0 * mu[r - n] + rs * z[r]) / (2.0 + rs);
            project_simplex(w);
            for (int r = seg_start[s]; r < seg_start[s + 1]; ++r) yg[r] = w[r - seg_start[s]];
        }
        for (int r = 0; r < rows; ++r) q[r] = rho[r] * (Ef[r] - yg[r] + lam[r]);
        ET.multiply(q, ETq);
        for (int v = 0; v < k; ++v) target[0][v] = f[v] - ETq[v] / sigma[v];
        tsddp.update_frequency_matrix(target);
        tsddp.solve();
        f = tsddp.frequencies[0];
        E.multiply(f, Ef);
        for (int r = 0; r < rows; ++r) lam[r] += Ef[r] - yg[r];
        ++t;

        // The certificate costs about as much as an iteration; evaluate it every check_every iterations.
        if (t % std::max(1, opt.check_every) != 0 && t < opt.max_iterations) continue;
        for (int r = 0; r < rows; ++r) beta[r] = rho[r] * lam[r];
        fd = dilute(f);
        L = loss(fd);
        d = dual_value(beta);
        if (opt.verbose) spdlog::info("iteration {}: loss {}, lower bound {}, gap {}", t, L, d, L - d);
        // L(u) is +inf while an SNV with v_i > 0 has y_i(u) = 0; never stop then
        if (std::isfinite(L) && L - d <= opt.eps) { res.converged = true; break; }
    }

    res.frequencies = fd;
    res.usages = u_from_f(fd);
    res.loss = L;
    res.lower_bound = d;
    res.iterations = t;
    res.runtime_ms = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start).count();
    return res;
}
};
