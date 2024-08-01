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

    helper.resize((n+2)*(K+2));
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

    Solver solver(n+1,K);
    solver.init(var,ref,link_list,r);

    real answer = solver.main(0.8,1e-6);

    printf("%.12lf\n",answer);

    for (int i = 0; i < n; i++){
        printf("%d : %.12lf\n",i, solver.F[i]);
    }

    double obj_recal = 0;
    printf("%.12lf\n", obj_recal);
    return 0;
}
