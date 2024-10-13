import sys
import json
import time
import gurobipy as gp
from gurobipy import GRB
import numpy as np
from math import log

def log_eps(x,eps=1e-6,s_n=3):
    if x < eps:
        return_val = log(eps)
        iter_val = (eps-x)/eps
        for i in range(1,s_n):
            return_val -= (iter_val**i)/i
        return return_val
    return log(x)

def mll(V,R,T,env,pass_t):
    tot_answer = 0.
    m_,n = V.shape
    result_F = np.zeros_like(V,dtype=float)
    n_intervals = pass_t
    m = 1
    M = np.zeros((1,n),dtype=float)
    
    model = gp.Model(env=env)

    f_vars = [[model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1) 
               for i in range(n)] 
              for p in range(m)]
    
    f_i_vars = [[[model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1) 
                  for break_point in range(n_intervals+1)] 
                 for i in range(n)]
                for p in range(m)]
    
    lin_expr = gp.LinExpr()
    lin_expr2 = gp.LinExpr()
    
    for p in range(m):
        for i in T[-1]:
            lin_expr += f_vars[p][i]
        model.addConstr(lin_expr<=1.0)
        lin_expr.clear()
    for p in range(m):
        for i in range(n):
            for j in T[i]:
                lin_expr += f_vars[p][j]
            model.addConstr(lin_expr<=f_vars[p][i])
            lin_expr.clear()
    
    update_constrs = [[None for i in range(n)] for p in range(m)]
    for p in range(m):
        for i in range(n):
            for k in range(n_intervals+1):
                lin_expr += f_i_vars[p][i][k];
                lin_expr2 += f_i_vars[p][i][k] # * break_points[p][i][k] placeholder only
            model.addConstr(lin_expr == 1)
            update_constrs[p][i] = model.addConstr(lin_expr2 == f_vars[p][i])
            lin_expr.clear()
            lin_expr2.clear()

    for _ in range(m_):
        M.fill(0.5)
        Range = 0.5
        for __ in range(1):
            L = M - Range
            U = M + Range
            break_points = [[[(1-break_point/n_intervals)*l+break_point/n_intervals*u 
                              for break_point in range(n_intervals+1) ]
                             for l,u in zip(ll,uu)]
                            for ll,uu in zip(L,U)]

            for p in range(m):
                for i in range(n):
                    f_vars[p][i].LB=L[p][i]
                    f_vars[p][i].UB=U[p][i]

            for p in range(m):
                for i in range(n):
                    for k in range(n_intervals+1):
                        model.chgCoeff(update_constrs[p][i], f_i_vars[p][i][k], newvalue=break_points[p][i][k])

            for p in range(m):
                for i in range(n):
                    for k in range(n_intervals+1):
                        lin_expr += f_i_vars[p][i][k] * (V[_][i]*log_eps(break_points[p][i][k])+
                                                  R[_][i]*log_eps(1-break_points[p][i][k]))
            
            model.setObjective(lin_expr, GRB.MAXIMIZE)
            
            lin_expr.clear()
            model.optimize()
            Range*=4/n_intervals
            for p in range(m):
                for i in range(n):
                    M[p][i] = f_vars[p][i].X
            for i in range(n):
                result_F[_][i] = M[0][i]
                    
        tot_answer += model.ObjVal

    return tot_answer,result_F

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.stderr.write("Usage: python " + sys.argv[0] + " <variant_matrix> <total_matrix> <tree> <nr_segments>\n")
        sys.exit(1)

    env = gp.Env(empty=True)
    env.setParam("OutputFlag",0)
    env.start()

    var = np.loadtxt(sys.argv[1])
    tot = np.loadtxt(sys.argv[2])
    ref = tot-var

    tree = []
    non_root = set([])
    with open(sys.argv[3],"r") as fin:
        for line in fin:
            if line[0] == "#":
                continue
            tree.append([int(_) for _ in line.split()[1:]])
            non_root.update(tree[-1])
    root, = list(set(range(len(tree))).difference(non_root))
    tree.append([root])

    time_start = time.time()


    obj,F = mll(var,ref,tree,env,int(sys.argv[4]))
    runtime = (time.time()-time_start)*1000
            
    results = {"objective":-obj, "frequency_matrix":F.tolist(), "runtime":runtime}
    print(json.dumps(results))
