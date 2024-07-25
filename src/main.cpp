#include <fstream>
#include <iostream>

#include "Solver.h"


int main(int argc, char ** argv) {
    std::vector<int> var, ref, dep;
    std::vector<std::list<int> > link_list;
    int n,K,r;

    std::ifstream fin(argv[1]);
    K = std::stoi(argv[2]);
    fin >> n;
    var.resize(n);
    ref.resize(n);
    dep.resize(n);
    link_list.resize(n);
    for (int i = 0; i < n; i++){
        fin >> var[i] ;
    }
    for (int i = 0; i < n; i++){
        fin >>  dep[i];
        ref[i] = dep[i] - var[i];
    }
    fin >> r;
    int u,v;
    for (int i = 1; i < n; i++){
        fin >> u >> v;
        link_list[u].push_back(v);
    }

    Solver solver(var,ref,K,link_list,r);
//    real answer = solver.main();

    solver.init_range(0,1);

    solver.dfs(solver.root);

    real answer = solver.answer();

    solver.backtrace();

    printf("%.12lf\n",answer);
    return 0;
}
