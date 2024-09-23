import json
import time
import argparse
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from gurobipy import *

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--clonal-matrix', type=str, required=False,
                        help='Numpy TXT file containing the input matrix B.')
    parser.add_argument('--tree', type=str, required=False,
                        help='Adjacency list describing the input tree.')
    parser.add_argument('--frequency-matrix', type=str, required=True,
                        help='Numpy TXT file containing the input frequency matrix F.')
    parser.add_argument('--weight-matrix', type=str, required=False,
                        help='Numpy TXT file containing the input weight matrix W.')
    parser.add_argument('--output', type=str, required=False,
                        help="Prefix of output files.")

    return parser.parse_args()

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

def one_vafpp_linear_program(B, F, W=None):
    if W is None:
        W = np.ones_like(F)

    if W.ndim == 1:
        W = W.reshape(1, -1)

    if F.ndim == 1:
        F = F.reshape(1, -1)

    # m = number of samples, n = number of clones
    n, m = F.shape[1], F.shape[0]

    model = Model("Primal QP")

    U = model.addMVar(shape=(m, n), lb=0, vtype=GRB.CONTINUOUS, name="U")
    Z = model.addMVar(shape=(m, n), lb=float('-inf'), vtype=GRB.CONTINUOUS, name="Z")

    model.addConstr(Z == F - U @ B)
    #model.addConstr(Z >= U @ B - F)

    for k in range(m):
        model.addConstr(U[k, :].sum() <= 1)

    model.setObjective(quicksum(Z[k, i] * Z[k, i] * W[k, i] for k in range(m) for i in range(n)), GRB.MINIMIZE)
    model.optimize()

    # grab U from the model
    U_hat = np.array([U[k, :].X for k in range(m)])

    return model.objVal, model.Runtime, U_hat

if __name__ == '__main__':
    args = parse_args()

    if args.tree:
        tree = nx.read_adjlist(args.tree, nodetype=int, create_using=nx.DiGraph())
        B = construct_clonal_matrix(tree)
    else:
        B = np.loadtxt(args.clonal_matrix)

    F = np.loadtxt(args.frequency_matrix)

    if args.weight_matrix:
        W = np.loadtxt(args.weight_matrix)
    else:
        W = None

    start = time.time()
    obj1, gurobi_time, U_hat = one_vafpp_linear_program(B, F, W=W)
    end = time.time()

    # print(U_hat) print these more nicely
    print(U_hat @ B)



