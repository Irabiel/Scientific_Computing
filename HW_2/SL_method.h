//
// Created by Matt Blomquist on 9/29/22.
//

#ifndef F22_MATH_233_LABS_SL_METHOD_H
#define F22_MATH_233_LABS_SL_METHOD_H

#include "Grid2d.h"
#include "math_tools.h"
#include <vector>

// Semi-Langrangian Method
class SL_method {
private:
    Grid2d sl_grid;
    std::vector<double> sol;
    std::vector<double> sol_old;
    std::vector<double> sol_true;
    std::vector<double> vel_u;
    std::vector<double> vel_v;
    std::vector<double> temp_vel;
    int nx;
    int ny;
    void find_trajectory(int n, double & x_d, double & y_d, double dt);
    void find_trajectoryRK2(int n, double & x_d, double & y_d, double dt);
    void sol_IC(std::vector<double> & sol0);
    void compute_vel(double x, double y);

public:
    SL_method(); // constructor
    void set_grid(Grid2d & new_grid); // set grid
    std::vector<double> get_sol(){ return sol; }        // access solution
    std::vector<double> get_True_sol(){ return sol_true; }        // access true solution
    void set_velocity(std::vector<double> & vel_u0, std::vector<double> & vel_v0);
    void set_velocity();
    void Solver(double dt, double Tf);
    double compute_error();

};


#endif //F22_MATH_233_LABS_SL_METHOD_H
