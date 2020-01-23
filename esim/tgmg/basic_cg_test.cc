#include "poisson_1d.hh"

int main() {
    poisson_1d po(11);
    po.init();
    po.solve(true);
    po.print_solution();
}
