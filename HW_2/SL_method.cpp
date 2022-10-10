//
// Created by Matt Blomquist on 9/29/22.
//

#include "SL_method.h"



void SL_method::find_trajectory(int n, double &x_d, double &y_d, double dt) {
    // RK1 Euler Method
    double x_0 = sl_grid.x_from_n(n);
    double y_0 = sl_grid.y_from_n(n);
    x_d = x_0 - dt * vel_u[n];
    y_d = y_0 - dt * vel_v[n];
}

void SL_method::find_trajectoryRK2(int n, double &x_d, double &y_d, double dt) {
    // RK1 Euler Method
    double x_0 = sl_grid.x_from_n(n);
    double y_0 = sl_grid.y_from_n(n);
    double x_temp = 0.;
    double y_temp = 0.;
    x_d = x_0 - dt * vel_u[n];
    y_d = y_0 - dt * vel_v[n];
}