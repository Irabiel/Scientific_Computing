//
// Created by Matt Blomquist on 9/29/22.
//

#ifndef F22_MATH_233_LABS_SL_METHOD_H
#define F22_MATH_233_LABS_SL_METHOD_H

#include "Grid2d.h"
#include <vector>

// Semi-Langrangian Method
class SL_method {
private:
    Grid2d sl_grid;
    std::vector<double> sol;
    std::vector<double> vel_u;
    std::vector<double> vel_v;
    void find_trajectory(int n, double & x_d, double & y_d, double dt);
    void find_trajectoryRK2(int n, double & x_d, double & y_d, double dt);
    void reset_temp();

public:
    void set_grid(Grid2d & new_grid){sl_grid = new_grid;} // set grid
    std::vector<double> get_sol(){ return sol; }        // access solution
    void set_velocity(std::vector<double> & vel_u0, std::vector<double> & vel_v0); // for constant velocity

};


#endif //F22_MATH_233_LABS_SL_METHOD_H
