//
// Created by irabiel on 10/9/22.
//

#include <iostream>
#include "math_tools.h"
#include "Grid2d.h"
#include <vector>
#include <cmath>
#include "FullMatrix.h"

using namespace std;

int main(){
    long N = 3;
    long M = 3;

    double xmin = -1.;
    double xmax = 1.;
    double ymin = -1.;
    double ymax = 1.;

    Grid2d newGrid(N, M, xmin, xmax, ymin, ymax);

   cout << minmod(-3.2 , -3.) << endl;

}