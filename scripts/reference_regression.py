import argparse
import json
import sys
import time
import numpy as np
import networkx as nx
import cvxpy as cp

def solve_cvxpy(V, R, tree, solver="ECOS", eps=1e-7, loss="binomial"):
    if len(V.shape) == 1:
        V = V.reshape(1, -1)
        print(V)
    m, n = V.shape
    f = cp.Variable((m, n))
    if loss == "l2":
        constraints = [f >= 0, f <= 1]
    else:
        constraints = [f >= eps, f <= 1 - eps]
    children = [[] for _ in range(n)]
    for node in tree.nodes():
        children[node] = list(tree.successors(node))
    for j in range(n):
        if children[j]:
            constraints.append(cp.sum(f[:, children[j]], axis=1) <= f[:, j])
    root_nodes = [node for node in tree.nodes() if tree.in_degree(node) == 0]
    for r_ in root_nodes:
        constraints.append(f[:, r_] <= 1.0)
    T = V + R
    if loss == "l2":
        F = np.divide(V, T, out=np.zeros_like(V), where=T != 0)
        obj = cp.sum_squares(f - F)
    else:
        obj = -cp.sum(cp.multiply(V, cp.log(f)) + cp.multiply(R, cp.log(1 - f)))
    prob = cp.Problem(cp.Minimize(obj), constraints)
    prob.solve(solver=solver, verbose=True)
    return {
        "status": prob.status,
        "objective": prob.value,
        "runtime": prob.solver_stats.solve_time # don't measure compilation time
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("variant_matrix")
    parser.add_argument("total_matrix")
    parser.add_argument("tree")
    parser.add_argument("-s", "--solver", choices=["ECOS", "CLARABEL", "MOSEK", "GUROBI"], default="ECOS")
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
