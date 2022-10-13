//
// Created by Matt Blomquist on 9/29/22.
//

#include "SL_method.h"
#include "cmath"

using namespace std;

SL_method::SL_method(){

}

void SL_method::compute_vel(double x, double y){
    temp_vel.assign(2, 0.);
    temp_vel[0] = -y;
    temp_vel[1] = x;
}

void SL_method::sol_IC(){
    sol_old.assign(nx*ny, 0.);

    for (int i = 0; i < nx*ny; i++){
        sol_old[i] = sqrt((sl_grid.x_from_n(i) - 0.25)*(sl_grid.x_from_n(i) - 0.25) + sl_grid.y_from_n(i)*sl_grid.y_from_n(i)) - 0.2;
    }
}

void SL_method::sol_True(){
    sol_true.assign(nx*ny, 0.);

    for (int i = 0; i < nx*ny; i++){
        sol_true[i] = sqrt((sl_grid.x_from_n(i) - 0.25)*(sl_grid.x_from_n(i) - 0.25) + sl_grid.y_from_n(i)*sl_grid.y_from_n(i)) - 0.2;
    }
}

void SL_method::update_sol_old(){
    for (int i = 0; i < nx*ny ; i++){
        sol_old[i] = sol[i];
    }
}

void SL_method::set_velocity() {
    for (int i = 0; i< nx* ny; i++){
        vel_u[i] = - sl_grid.y_from_n(i) ;
        vel_v[i] = sl_grid.x_from_n(i) ;
    }
}

void SL_method::find_trajectory(int n, double &x_d, double &y_d, double dt) {
    // RK1 Euler Method
    double x_0 = sl_grid.x_from_n(n);
    double y_0 = sl_grid.y_from_n(n);
    x_d = x_0 - dt * vel_u[n];
    y_d = y_0 - dt * vel_v[n];
}

void SL_method::find_trajectoryRK2(int n, double &x_d, double &y_d, double dt) {
    // RK2 Euler Method
    double x_0 = sl_grid.x_from_n(n);
    double y_0 = sl_grid.y_from_n(n);
    compute_vel(x_0, y_0);
    x_star = x_0 - .5*dt * temp_vel[0];
    y_star = y_0 - .5*dt * temp_vel[1];
    compute_vel(x_star, y_star);
    x_d = x_0 - dt * temp_vel[0];
    y_d = y_0 - dt * temp_vel[1];
}

void SL_method::Solver(double dt, double Tf) {

    sol.assign(nx * ny, 0.);
    sol_old.assign(nx * ny, 0.);

    std::vector<double> xd;
    std::vector<double> yd;
    xd.assign(nx, 0.);
    yd.assign(ny, 0.);

    set_velocity();
    for (int i = 0; i < nx ; i++)
        for (int j = 0; j < ny; j++)
            find_trajectory(i, xd[i], yd[j], dt);

    int nt = ceil(Tf / dt);
    sol_IC();

    for (int n = 0; n < nt; n++) {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                sol[i*nx + j] = quadratic_interpolation(sl_grid, sol_old, xd[i], yd[j]);
        update_sol_old();
    }

}

double SL_method::compute_error(){
    double err = 0.;
    sol_True();

    for (int i = 0; i < nx*ny; i++)
        err += abs(sol[i] - sol_true[i]);

    return err;
}

void SL_method::set_grid(Grid2d & new_grid){
    sl_grid = new_grid;
    nx = sl_grid.get_numb_x();
    ny = sl_grid.get_numb_y();
} // set grid