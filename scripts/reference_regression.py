import argparse
import json
import sys
import time
import numpy as np
import networkx as nx
import cvxpy as cp

def solve_cvxpy(V, R, tree, solver="ECOS", eps=1e-4, loss="binomial"):
    if len(V.shape) == 1:
        V = V.reshape(1, -1)

    m, n = V.shape

    f = cp.Variable(n)
    u = cp.Variable(n)

    B = construct_clonal_matrix(tree)
    if loss == "l2":
        constraints = [u >= 0, cp.sum(u) <= 1, cp.matmul(u, B) == f]
    else:
        constraints = [u >= 0, cp.sum(u) <= 1, cp.matmul(u, B) == f, f >= eps, f <= 1 - eps]

    obj_value = 0
    runtime = 0
    failed_subproblem = False
    if loss == "l2":
        F_obs = cp.Parameter(n)
        obj = cp.sum_squares(f - F_obs)
        prob = cp.Problem(cp.Minimize(obj), constraints)

        for i in range(m):
            T = V[i,:] + R[i,:]
            F = np.divide(V[i], T, out=np.zeros_like(V[i]), where=T != 0)
            F_obs.value = F
            try:
                prob.solve(solver=solver)
            except cp.error.SolverError:
                failed_subproblem = True
                break
            obj_value += prob.value
            runtime += prob.solver_stats.solve_time
    else:
        V_obs = cp.Parameter(n, nonneg=True)
        R_obs = cp.Parameter(n, nonneg=True)
        obj = -cp.sum(cp.multiply(V_obs, cp.log(f)) + cp.multiply(R_obs, cp.log(1 - f)))
        prob = cp.Problem(cp.Minimize(obj), constraints)

        for i in range(m):
            V_obs.value = V[i]
            R_obs.value = R[i]
            try:
                if solver == 'CLARABEL':
                    prob.solve(solver=solver, verbose=False, max_iter=200)
                elif solver == 'MOSEK':
                    prob.solve(solver=solver, verbose=False, accept_unknown=True)
                elif solver == 'ECOS':
                    prob.solve(solver=solver, verbose=False, max_iters=400)
                elif solver == 'SCS':
                    prob.solve(solver=solver, verbose=True)

            except cp.error.SolverError:
                failed_subproblem = True
                break
            obj_value += prob.value
            runtime += prob.solver_stats.solve_time

    return {
        "failed_subproblem": failed_subproblem,
        "objective": obj_value,
        "runtime": runtime * 1000 # convert to milliseconds
    }

"""
Constructs the clonal matrix from a clonal tree.
"""
def construct_clonal_matrix(tree):
    n = len(tree.nodes)

    B = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if nx.has_path(tree, i, j):
                B[j, i] = 1

    return B

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("variant_matrix")
    parser.add_argument("total_matrix")
    parser.add_argument("tree")
    parser.add_argument("-s", "--solver", choices=["ECOS", "CLARABEL", "MOSEK", "GUROBI", "SCS"], default="ECOS")
    parser.add_argument("-l", "--loss", choices=["binomial", "l2"], default="binomial")
    args = parser.parse_args()

    V = np.loadtxt(args.variant_matrix)
    T = np.loadtxt(args.total_matrix)
    V = V
    R = T - V

    tree = nx.read_adjlist(args.tree, nodetype=int, create_using=nx.DiGraph())
    result = solve_cvxpy(V, R, tree, args.solver, 1e-7, args.loss)
    print(json.dumps(result))

if __name__ == "__main__":
    main()
