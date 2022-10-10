//
// Created by Matt Blomquist on 9/27/22.
//

#include "math_tools.h"
#include <cmath>

// suppose we have grid in 2D [xmin,xmax] x [ymin,ymax]
// find cell in which (x,y) belongs
// find weighted avg of values (?)

double bilinear_interpolation(Grid2d & grid,std::vector<double> & func,double x, double y){
    //  get xmin and get ymin
    // if (x,y) are outside domain, throw warning/error message
    if (x < grid.get_xmin() || x > grid.get_xmax() || y < grid.get_ymin() || y > grid.get_ymax() ) {
        throw std::invalid_argument("ERROR: (x,y) is outside of grid");
    }

    double phi;
    double dx = grid.get_dx();
    double dy = grid.get_dy();

    std::cout <<"dx: " << dx << " dy: " << dy << std::endl;

    // get i and j to find which cell C x belongs to
    //  i = floor((x - xmin)/dx)
    //  j = floor((y-ymin)/dx)
    int i = floor( (x - grid.get_xmin()) / dx);
    int j = floor( (y - grid.get_ymin()) / dy);

    std::cout << "i: " << i << " j: " << j << std::endl;

    double x_i = grid.get_xmin() + i * dx;
    double y_j = grid.get_ymin() + j * dy;
    double x_ip1 = x_i + dx;
    double y_jp1 = y_j + dy;

    // Use quadratic interpolation (formula in Lab 3) to get value at x
    // (i.e. think weighted avg)
    // (i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1) are the corners of the cell C
    phi  = func[grid.n_from_ij(i    ,   j  )]  * ( x_ip1 - x   ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1,   j  )]  * ( x     - x_i ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i    , j+1)]  * ( x_ip1 - x   ) * ( y     - y_j ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1, j+1)]  * ( x     - x_i ) * ( y     - y_j ) / (dx*dy) ;

    return phi;
}

double quadratic_interpolation(Grid2d & grid,std::vector<double> & func,double x, double y){
    if (x < grid.get_xmin() || x > grid.get_xmax() || y < grid.get_ymin() || y > grid.get_ymax() ) {
        throw std::invalid_argument("ERROR: (x,y) is outside of grid");
    }

    double phi;
    double dx = grid.get_dx();
    double dy = grid.get_dy();

    double ax;
    double ay;
    double bx;
    double by;

    std::cout <<"dx: " << dx << " dy: " << dy << std::endl;

    // get i and j to find which cell C x belongs to
    int i = floor( (x - grid.get_xmin()) / dx);
    int j = floor( (y - grid.get_ymin()) / dy);

    // check if the x component is on either boundary or next to the left boundary
    if (i == 0) {
        FDd2(func[grid.n_from_ij(i + 1, j)], func[grid.n_from_ij(i, j)], func[grid.n_from_ij(i - 1, j)], dx);
        ax = FDd2(func[grid.n_from_ij(i + 2, j)], func[grid.n_from_ij(i + 1, j)], func[grid.n_from_ij(i, j)],
                         dx);
        bx = FDd2(func[grid.n_from_ij(i + 2, j)], func[grid.n_from_ij(i + 1, j)], func[grid.n_from_ij(i, j)],
                         dx);
    }
    else if (i == grid.get_numb_x() -1){
        ax = FDd2(func[grid.n_from_ij(i, j)], func[grid.n_from_ij(i - 1, j)], func[grid.n_from_ij(i-2, j)],
                         dx);
        bx = FDd2(func[grid.n_from_ij(i, j)], func[grid.n_from_ij(i - 1, j)], func[grid.n_from_ij(i-2, j)],
                         dx);
    }
    else if (i == grid.get_numb_x() -2){
        ax = FDd2(func[grid.n_from_ij(i + 1, j)], func[grid.n_from_ij(i, j)], func[grid.n_from_ij(i - 1, j)], dx);
        bx = FDd2(func[grid.n_from_ij(i + 1, j)], func[grid.n_from_ij(i, j)], func[grid.n_from_ij(i - 1, j)], dx);
    }
    else{
        ax = FDd2(func[grid.n_from_ij(i + 1, j)], func[grid.n_from_ij(i, j)], func[grid.n_from_ij(i - 1, j)], dx);
        bx = FDd2(func[grid.n_from_ij(i + 2, j)], func[grid.n_from_ij(i + 1, j)], func[grid.n_from_ij(i, j)], dx);
    }

    // check if the y component is on either boundary or next to the left boundary
    if (j == 0) {
        ax = FDd2(func[grid.n_from_ij(i, j + 2)], func[grid.n_from_ij(i, j + 1)], func[grid.n_from_ij(i, j)],
                         dy);
        bx = FDd2(func[grid.n_from_ij(i, j + 2)], func[grid.n_from_ij(i, j + 1)], func[grid.n_from_ij(i, j)],
                         dy);
    }
    else if (j == grid.get_numb_y()-1){
        ax = FDd2(func[grid.n_from_ij(i, j)], func[grid.n_from_ij(i, j -1)], func[grid.n_from_ij(i, j - 2)],
                         dy);
        bx = FDd2(func[grid.n_from_ij(i, j)], func[grid.n_from_ij(i, j - 1)], func[grid.n_from_ij(i, j - 2)],
                         dy);
    }
    else if (j == grid.get_numb_y()-2){
        ax = FDd2(func[grid.n_from_ij(i, j + 1)], func[grid.n_from_ij(i, j)], func[grid.n_from_ij(i, j - 1)], dy);
        bx = FDd2(func[grid.n_from_ij(i, j + 1)], func[grid.n_from_ij(i, j)], func[grid.n_from_ij(i, j - 1)], dy);
    }
    else{
        ax = FDd2(func[grid.n_from_ij(i, j + 1)], func[grid.n_from_ij(i, j)], func[grid.n_from_ij(i, j - 1)], dy);
        bx = FDd2(func[grid.n_from_ij(i, j + 2)], func[grid.n_from_ij(i, j + 1)], func[grid.n_from_ij(i, j)], dy);
    }

    std::cout << "i: " << i << " j: " << j << std::endl;

    double x_i = grid.get_xmin() + i * dx;
    double y_j = grid.get_ymin() + j * dy;
    double x_ip1 = x_i + dx;
    double y_jp1 = y_j + dy;

    // Use quadratic interpolation (formula in Assignment 3) to get value at x
    // (i.e. think weighted avg)
    // (i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1) are the corners of the cell C
    phi  = func[grid.n_from_ij(i    ,   j  )]  * ( x_ip1 - x   ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1,   j  )]  * ( x     - x_i ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i    , j+1)]  * ( x_ip1 - x   ) * ( y     - y_j ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1, j+1)]  * ( x     - x_i ) * ( y     - y_j ) / (dx*dy) ;
    phi -= ( x     - x_i ) * ( x_ip1 - x   ) * minmod(ax, bx) ;
    phi -= ( y     - y_j ) * ( y_jp1 - y   ) * minmod(ay, by) ;
    return phi;
}

inline double FDd2(double fxp,double fx,double fxm, double dx){
    return ( (fxp - 2*fx + fxm) / (dx*dx) );
}

double minmod(double a, double b){
    if (a*b < 0)
        return 0;
    else if ( abs(a) < abs(b))
        return a;

    return b;

}


