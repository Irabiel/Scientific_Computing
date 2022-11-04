//
// Created by Matt Blomquist on 9/29/22.
//

#include "SL_method.h"
#include "cmath"
#include <random>

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

    std::vector<double> xd, yd;
    xd.assign(nx*ny, 0.);
    yd.assign(ny*nx, 0.);

    for (int i = 0; i < nx ; i++)
        for (int j = 0; j < ny; j++)
            find_trajectoryRK2(i + j*ny, xd[i + j*ny], yd[i + j*ny], dt);

    int nt = floor(Tf / dt);

    sol.assign(nx * ny, 0.);
    sol_IC(sol_old);

    char file_name[250];
    sprintf(file_name,"../q3VTK.vtk");
    sl_grid.print_VTK_format(file_name);

    for (int n = 0; n < nt; n++) {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                sol[i + j*ny] = quadratic_interpolation(sl_grid, sol_old, xd[i + j*ny], yd[i + j*ny]);

        if (REI == 1)
            Reinitialization_Equation( sl_grid.get_dx()/10, sol);

        sol_old = sol;

        char data_name[250];
        if (n % 10 == 0 || n == nt-1){
            sprintf(data_name,"num_iter=%d", n);
            sl_grid.print_VTK_format(sol, data_name,file_name);
        }

    }
    sol_IC(sol_true);
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

void SL_method::Godunov_Scheme(std::vector<double> & G) {
    G.assign(nx*ny, 0.);
    double a, b, c, d, ap, bp, cp, dp, am, bm, cm, dm;
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            a = (i == 0   ) ? 0.0    : (RIE_sol[sl_grid.n_from_ij(i,      j)] - RIE_sol[sl_grid.n_from_ij(i-1, j)]  )/sl_grid.get_dx();
            b = (i == nx-1) ? 0.0    : (RIE_sol[sl_grid.n_from_ij(i+1, j)] - RIE_sol[sl_grid.n_from_ij(i,      j)]  )/sl_grid.get_dx();
            c = (j == 0   ) ? 0.0    : (RIE_sol[sl_grid.n_from_ij(i,      j)] - RIE_sol[sl_grid.n_from_ij(i,   j-1)])/sl_grid.get_dy();
            d = (j == ny-1) ? 0.0    : (RIE_sol[sl_grid.n_from_ij(i, j+1)] - RIE_sol[sl_grid.n_from_ij(i,      j)]  )/sl_grid.get_dy();

            if (sol[sl_grid.n_from_ij(i, j)] > 0.){
                ap = max(a, 0.);
                bm = min(b, 0.);
                cp = max(c, 0.);
                dm = min(d, 0.);
                G[sl_grid.n_from_ij(i, j)] = sqrt( max(ap*ap, bm*bm) + max(cp*cp, dm*dm) ) - 1.;
            }
            else if (sol[sl_grid.n_from_ij(i, j)] < 0.){
                am = min(a, 0.);
                bp = max(b, 0.);
                cm = min(c, 0.);
                dp = max(d, 0.);
                G[sl_grid.n_from_ij(i, j)] = sqrt( max(am*am, bp*bp) + max(cm*cm, dp*dp) ) - 1.;
            }
            else
                G[sl_grid.n_from_ij(i, j)] = 0.0;
        }

    }
}

void SL_method::Reinitialization_Equation(double dt, std::vector<double> & Sol) {
    double epsilon = max(sl_grid.get_dx(), sl_grid.get_dy());
    S.assign(nx*ny, 0.);
    G.assign(nx*ny, 0.);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            S[sl_grid.n_from_ij(i, j)] = (Sol[sl_grid.n_from_ij(i, j)]
                                          / sqrt(Sol[sl_grid.n_from_ij(i, j)] * Sol[sl_grid.n_from_ij(i, j)] +
                                                 epsilon * epsilon));
        }
    }
    RIE_sol = Sol;

    for (int k = 0; k < 10; k++){
        Godunov_Scheme(G);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                RIE_sol[sl_grid.n_from_ij(i, j)] -= dt * S[sl_grid.n_from_ij(i, j)] * G[sl_grid.n_from_ij(i, j)];
            }
        }
    }
    Sol = RIE_sol;
}

void SL_method::Reinitialization_Equation_with_plot(double dt, std::vector<double> & Sol, int num_iter) {
    double epsilon = max(sl_grid.get_dx(), sl_grid.get_dy());
    S.assign(nx*ny, 0.);
    G.assign(nx*ny, 0.);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            S[sl_grid.n_from_ij(i, j)] = (Sol[sl_grid.n_from_ij(i, j)]
                                          / sqrt(Sol[sl_grid.n_from_ij(i, j)] * Sol[sl_grid.n_from_ij(i, j)] +
                                                 epsilon * epsilon));
        }
    }
    RIE_sol = Sol;

    char file_name[250];
    sprintf(file_name,"../q2VTK.vtk");
    sl_grid.print_VTK_format(file_name);

    char data_name[250];
    sprintf(data_name,"num_iter=%d", 0);
    sl_grid.print_VTK_format(RIE_sol, data_name,file_name);

    for (int k = 0; k < num_iter; k++){
        Godunov_Scheme(G);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                RIE_sol[sl_grid.n_from_ij(i, j)] -= dt * S[sl_grid.n_from_ij(i, j)] * G[sl_grid.n_from_ij(i, j)];
            }
        }

        char data_name[250];
        sprintf(data_name,"num_iter=%d", k+1);
        sl_grid.print_VTK_format(RIE_sol, data_name,file_name);

    }
    Sol = RIE_sol;
}

void SL_method::error_map(vector<double> &sol, vector<double> &True, vector<double> &error) {
    error.assign(nx*ny,0.);
    for (int n = 0; n < nx*ny; n++) {
        error[n] = sol[n] - True[n];
    }
}

void SL_method::Add_nosie(std::vector<double> & vec){

    for (int i = 0; i < nx*ny; ++i)
        if (abs(vec[i]) > 0.02)
            vec[i] +=  .01 * ((float)rand() / (float)RAND_MAX);

}