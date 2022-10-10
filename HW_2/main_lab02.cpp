#include <iostream>
#include <vector>
#include "Grid2d.h"
#include "FullMatrix.h"



// using namespace std;  < -- Streamline?

int main() {

    std::cout << "Hello, Lab 02" << std::endl;

    FullMatrix A(5);
    //    A.add_Element(6,1,3);     // check throws for invalid size
    A.add_Element(0, 0, 2);
    A.add_Element(1, 1, 3);
    A.add_Element(2, 3, 5);

    std::cout << "A Matrix:" << std::endl;
    A.display();

    std::vector<double> x, y;

//    x.resize(5);
//    y.resize(5);
//
//    x.assign(5, 1.);

//    y = A.mat_Vec_Prod(x);

//    for (int i = 0; i < 5; i++)
//        std::cout << y[i] << std::endl;

    // Create an Identity Matrix
    FullMatrix B(5);
    B.add_Element(0,0,1);
    B.add_Element(1,1,1);
    B.add_Element(2,2,1);
    B.add_Element(3,3,1);
    B.add_Element(4,4,1);

    FullMatrix AB(5);
    A.mat_Mat_Prod(B, AB);

    std::cout << "AB Matrix:" << std::endl;
    AB.display();


    return 0;
} 