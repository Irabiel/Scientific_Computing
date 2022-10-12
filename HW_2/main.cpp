//
// Created by irabiel on 10/9/22.
//

#include <iostream>
#include "math_tools.h"
#include "Grid2d.h"
#include <vector>
#include <cmath>
#include "FullMatrix.h"
#include "SL_method.h"

using namespace std;

int main(){
    long N = 3;
    long M = 3;

    double xmin = -1.;
    double xmax = 1.;
    double ymin = -1.;
    double ymax = 1.;

    Grid2d newGrid(N, M, xmin, xmax, ymin, ymax);

    SL_method SLM = SL_method();

    SLM.set_grid(newGrid);

    SLM.Solver(.1,1);

    std:: cout << SLM.compute_error() << std::endl;


}