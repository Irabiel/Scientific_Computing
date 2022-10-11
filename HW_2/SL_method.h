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
    std::vector<double> vel_u;
    std::vector<double> vel_v;
    std::vector<double> temp_vel;
    double x_star;
    double y_star;
    void find_trajectory(int n, double & x_d, double & y_d, double dt);
    void find_trajectoryRK2(int n, double & x_d, double & y_d, double dt);
    void sol_IC();
    void update_sol_old();
    std::vector<double> compute_vel(double x, double y);

public:
    void set_grid(Grid2d & new_grid){sl_grid = new_grid;} // set grid
    std::vector<double> get_sol(){ return sol; }        // access solution
    void set_velocity(std::vector<double> & vel_u0, std::vector<double> & vel_v0); // for constant velocity
    void Solver(double dt, double Tf);
};


#endif //F22_MATH_233_LABS_SL_METHOD_H
