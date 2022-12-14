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
    long M = N;

    double xmin = -1.;
    double xmax = 1.;
    double ymin = -1.;
    double ymax = 1.;

    Grid2d newGrid(N, M, xmin, xmax, ymin, ymax);

    double dx = newGrid.get_dx();
    double dy = newGrid.get_dy();

    double ratio = .5;
    double dt = max(dx, dy) / ratio;

    SL_method SLM = SL_method();
    SLM.set_grid(newGrid);
    SLM.set_REI(1);
    SLM.Solver(dt, 2 * M_PI);

    vector<double> sol = SLM.get_sol();
    SLM.get_True_sol();
    vector<double> True_sol = SLM.get_True_sol();

    newGrid.print_VTK_format("../SLM_approx.vtk");
    newGrid.print_VTK_format(sol, "Approx", "../SLM_approx.vtk");
    newGrid.print_VTK_format("../SLM_True.vtk");
    newGrid.print_VTK_format(True_sol, "Approx", "../SLM_True.vtk");

    vector<double> err_map;
    SLM.error_map(sol, True_sol,err_map);
    newGrid.print_VTK_format("../Error_map.vtk");
    newGrid.print_VTK_format(err_map, "Error_map", "../Error_map.vtk");

    SLM.Add_nosie(True_sol);
    SLM.Reinitialization_Equation_with_plot( newGrid.get_dx()/10, True_sol, 100);

    return 0;
}