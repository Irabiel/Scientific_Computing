//
// Created by Matt Blomquist on 9/27/22.
//
#include <iostream>
#include "math_tools.h"
#include "Grid2d.h"
#include <vector>
#include <cmath>

int main() {
    std::cout << "Hello, Lab 3!  " << std::endl;
    int N = 40;
    double xmin = 0.;
    double xmax = 1.;
    double ymin = 0.;
    double ymax = 1.;
    Grid2d new_grid = Grid2d(N, N, xmin, xmax, ymin, ymax );

    std::vector<double> vals;
    // allocate space in memory for vals
    vals.assign(N*N, 1.);

    for (int i = 0; i < N*N; i++){
        // smooth function for interpolation (Scott's fav) lol
        vals[i] = cos(new_grid.x_from_n(i)) + cos(new_grid.y_from_n(i));
        //std::cout << vals[i] << std::endl;
    }



    // choose (x,y) to compare interpolation vs. truth
    double pos_x = 0.12225;
    double pos_y = 0.2225;
    double pos = cos(pos_x) + cos(pos_y);

    std::cout << "pos_x: " << pos_x << " pos_y: "  << pos_y << std::endl;
    std::cout << "pos: " << pos << std::endl;

    // print out error (interpolation - truth)
    std::cout << "Error: " << std::abs(bilinear_interpolation(new_grid, vals, pos_x, pos_y) - pos) << " " << minmod(-1,2) << std::endl;
    return 0;
}