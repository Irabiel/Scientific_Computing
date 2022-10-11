//
// Created by Matt Blomquist on 9/29/22.
//

#include "SL_method.h"
#include "cmath"


std::vector<double> compute_vel(double x, double y){
    std::vector<double> vel_temp ;
    vel_temp.assign(2, 0.);
    vel_temp[0] = -y;
    vel_temp[1] = x;
}

void SL_method::sol_IC(){
    int nx = sl_grid.get_numb_x();
    int ny = sl_grid.get_numb_y();
    sol_old.assign(nx*ny, 0.);

    for (int i = 0; i < (nx-1)*(ny-1); i++){
        sol_old[i] = sqrt((sl_grid.x_from_n(i) - 0.25)*(sl_grid.x_from_n(i) - 0.25) + sl_grid.x_from_n(i)*sl_grid.x_from_n(i)) - 0.2;
    }

}

void SL_method::update_sol_old(){
    for (int i = 0; i < sol_old.size() - 1; i++){
        sol_old[i] = sol[i];
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
    temp_vel = compute_vel(x_0, y_0);
    x_star = x_0 - .5*dt * temp_vel[0];
    y_star = y_0 - .5*dt * temp_vel[1];
    temp_vel = compute_vel(x_star, y_star);
    x_d = x_0 - dt * temp_vel[0];
    y_d = y_0 - dt * temp_vel[1];
}

void SL_method::Solver(double dt, double Tf) {

int nx = sl_grid.get_numb_x();
int ny = sl_grid.get_numb_y();
sol.assign(nx * ny, 0.);

std::vector<double> xd;
std::vector<double> yd;
xd.assign(nx, 0.);
yd.assign(ny, 0.);


for (int i = 0; i < (nx)*(ny); i++){
    find_trajectoryRK2( i, xd[i], yd[i],  dt);
}

int nt = ceil(Tf/dt);
sol_IC();

for (int n = 1; n < nt; n++){
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            sol[i+j] = quadratic_interpolation(sl_grid, sol_old, xd[i], yd[j]);
    update_sol_old();
}


}