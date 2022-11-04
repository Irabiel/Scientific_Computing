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
    std::vector<double> sol  , sol_old, sol_true;
    std::vector<double> vel_u, vel_v  , temp_vel;
    std::vector<double> G    , S      ,  RIE_sol;
    int nx, ny;
    int REI = 0;
    void find_trajectory(int n, double & x_d, double & y_d, double dt);
    void find_trajectoryRK2(int n, double & x_d, double & y_d, double dt);

    void compute_vel(double x, double y);

public:
    SL_method(); // constructor
    void sol_IC(std::vector<double> & sol0);
    void set_grid(Grid2d & new_grid); // set grid
    std::vector<double> get_sol(){ return sol; };     // access solution
    std::vector<double> get_True_sol(){ return sol_true; };      // access true solution
    void set_velocity(std::vector<double> & vel_u0, std::vector<double> & vel_v0);
    void Add_nosie(std::vector<double> & vec);
    void set_velocity();
    void set_REI(double reI){REI = reI;};
    void Solver(double dt, double Tf);
    double compute_error();
    void error_map(std::vector<double> & sol, std::vector<double> & True, std::vector<double> & error);
    void Reinitialization_Equation(double dt, std::vector<double> & Sol);
    void Reinitialization_Equation_with_plot(double dt, std::vector<double> & Sol, int num_iter);
    void Godunov_Scheme(std::vector<double> & G);
};




#endif //F22_MATH_233_LABS_SL_METHOD_H
