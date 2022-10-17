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

void SL_method::sol_IC(std::vector<double> & sol0){
    sol0.assign(nx*ny, 0.);

    for (int i = 0; i < nx*ny; i++){
        sol0[i] = sqrt((sl_grid.x_from_n(i) - 0.25)*(sl_grid.x_from_n(i) - 0.25) + sl_grid.y_from_n(i)*sl_grid.y_from_n(i)) - 0.2;
    }
}

void SL_method::set_velocity() {
    vel_v.assign(nx*ny,0.);
    vel_u.assign(nx*ny,0.);

    for (int i = 0; i< nx* ny; i++){
        vel_u[i] = - sl_grid.y_from_n(i) ;
        vel_v[i] = sl_grid.x_from_n(i) ;
    }
}

void SL_method::find_trajectory(int n, double &x_d, double &y_d, double dt) {
    // RK1 Euler Method
    double x_0 = sl_grid.x_from_n(n);
    double y_0 = sl_grid.y_from_n(n);
    compute_vel(x_0, y_0);
    x_d = x_0 - dt * temp_vel[0];
    y_d = y_0 - dt * temp_vel[1];
}

void SL_method::find_trajectoryRK2(int n, double &x_d, double &y_d, double dt) {
    // RK2 Euler Method
    double x_0 = sl_grid.x_from_n(n);
    double y_0 = sl_grid.y_from_n(n);
    compute_vel(x_0, y_0);
    double x_star = x_0 - .5*dt * temp_vel[0];
    double y_star = y_0 - .5*dt * temp_vel[1];
    compute_vel(x_star, y_star);
    x_d = x_0 - dt * temp_vel[0];
    y_d = y_0 - dt * temp_vel[1];
}

void SL_method::Solver(double dt, double Tf) {

    sol.assign(nx * ny, 0.);
    sol_old.assign(nx * ny, 0.);

    std::vector<double> xd;
    std::vector<double> yd;
    xd.assign(nx*ny, 0.);
    yd.assign(ny*nx, 0.);
    for (int i = 0; i < nx ; i++)
        for (int j = 0; j < ny; j++) {
            find_trajectoryRK2(i + j*ny, xd[i + j*ny], yd[i + j*ny], dt);
        }

    int nt = ceil(Tf / dt);
    cout << "number of time steps nt = " << nt << endl;

    sol_IC(sol_old);

    for (int n = 0; n < nt; n++) {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                sol[i + j*ny] = quadratic_interpolation(sl_grid, sol_old, xd[i + j*ny], yd[i + j*ny]);
        sol_old = sol;
    }
}

double SL_method::compute_error(){
    double err = 0.;
    sol_IC(sol_true);

    for (int i = 0; i < nx*ny; i++) {
        err = max(abs(sol[i] - sol_true[i]), err); // inf norm
//        err += abs(sol[i] - sol_true[i]); // one norm
    }

    return err;
}

void SL_method::set_grid(Grid2d & new_grid){
    sl_grid = new_grid;
    nx = sl_grid.get_numb_x();
    ny = sl_grid.get_numb_y();
} // set grid