#include <iostream>
#include <vector>
#include "Grid2d.h"

// using namespace std;  < -- Streamline?

int main() {
    std::cout << "Hello, MATH233! :)" << std::endl;

    std::cout << abs(-1.5) << std::endl;

    long N = 3;
    long M = 3;

    double xmin = -1.;
    double xmax = 1.;
    double ymin = 2.;
    double ymax = 6.;

    Grid2d newGrid(N, M, xmin, xmax, ymin, ymax);

    std::cout << "This is dx: " << newGrid.get_dx() << std::endl;
    std::cout << "This is dy: " << newGrid.get_dy() << std::endl;

    std::cout << "Test i_from_n(7): " << newGrid.i_from_n(7) << std::endl;
    std::cout << "Test j_from_n(7): " << newGrid.j_from_n(7) << std::endl;
    std::cout << "Test n_from_ij(2,2): " << newGrid.n_from_ij(2, 2) << std::endl;

    std::cout << "Test x_from_n(2): " << newGrid.x_from_n(2) << std::endl;
    std::cout << "Test y_from_n(2): " << newGrid.y_from_n(2) << std::endl;

    std::vector<double> values;
    values.resize(N*M);
    for (int i = 0; i < N*M; i++){
        values[i] = i;
    }

    newGrid.print_VTK_format("../newVTK.vtk");
    newGrid.print_VTK_format(values, "newData", "../newVTK.vtk");

    return 0;

}

