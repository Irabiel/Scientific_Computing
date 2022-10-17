//
// Created by irabiel on 10/9/22.
//

#include <iostream>
#include "math_tools.h"
#include "Grid2d.h"
#include <vector>
#include <cmath>
#include "SL_method.h"

using namespace std;

int main(){
    long N = 100;
    long M = 100;

    double xmin = -1.;
    double xmax = 1.;
    double ymin = -1.;
    double ymax = 1.;

    Grid2d newGrid(N, M, xmin, xmax, ymin, ymax);

    double dx = newGrid.get_dx();
    double dy = newGrid.get_dy();

    double dt = max(dx,dy) * 1.5;

    SL_method SLM = SL_method();

    SLM.set_grid(newGrid);

    SLM.Solver(dt,2*M_PI);

    double err = SLM.compute_error();

    std:: cout << err << std::endl;

    vector<double> sol = SLM.get_sol();
    SLM.get_True_sol();
    vector<double> True_sol = SLM.get_True_sol();

    newGrid.print_VTK_format("../SLM_approx.vtk");
    newGrid.print_VTK_format(sol, "Approx", "../SLM_approx.vtk");
    newGrid.print_VTK_format("../SLM_True.vtk");
    newGrid.print_VTK_format(True_sol, "Approx", "../SLM_True.vtk");

    return 0;
}