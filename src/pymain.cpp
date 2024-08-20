//
// Created by Yuanyuan Qi on 8/10/24.
//

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include <cassert>

#include "Solver.h"

class pymain {
public:
    Solver solver;

    double Target_F_eps, Frac;
    int N_intervals, init_n_intervals;

    pymain(int max_n, int k=10, real target_F_eps=1e-4, real frac=0.8);

    void set_n_intervals(int);
    void set_frac(real);
    void set_target(real);


    real optimize(std::vector<std::list<int> > &tree_in_lists_of_children, int root, std::vector<int> &var,
                  std::vector<int> &ref);

    const std::vector<real> & F() const;
};

namespace py=pybind11;

pymain::pymain(int max_n, int k, real target_F_eps, real frac):
        Target_F_eps(target_F_eps),
        N_intervals(k),
        init_n_intervals(k),
        Frac(frac),
        solver(max_n,k)
{
}

real pymain::optimize(std::vector<std::list<int> > & tree_in_lists_of_children, int root, std::vector<int> &var,
                        std::vector<int> &ref) {
    solver.init(var,ref,tree_in_lists_of_children,root);
    real answer = solver.main(Frac,Target_F_eps);
    return answer;
}

const std::vector<real> &pymain::F() const {
    return solver.F;
}

void pymain::set_n_intervals(int n_intervals) {
    assert(n_intervals <= init_n_intervals);
    N_intervals = n_intervals;
    solver.n_intervals = n_intervals;
}

void pymain::set_frac(real frac) {
    Frac = frac;
}

void pymain::set_target(real target) {
    Target_F_eps = target;
}

PYBIND11_MODULE(test, m) {
    py::class_<pymain>(m, "CVXTree")
            .def(py::init<int, int, real, real>(), py::arg("max_n"), py::arg("max_n_interval") = 10,
                 py::arg("F_eps") = 1e-4, py::arg("base") = 0.8)
            .def("optimize", &pymain::optimize, py::arg("Lists_of_children"), py::arg("root"),
                 py::arg("para1"), py::arg("para2"))
            .def("F", &pymain::F)
            .def("set_n_intervals", &pymain::set_n_intervals)
            .def("set_base", &pymain::set_frac)
            .def("set_F_eps",&pymain::set_frac);
}