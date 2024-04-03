#pragma once
#include <cmath>
#include <iostream>
#include <vector>
#include "Shape_functions.h"
#include "Shape_fct_1D.h"

class Volume_Matrix_Integral
{
    //
    // Handles integral over the volume dV with terms in 𝜑ᵢ×𝜑ⱼ
    // Used in the left-handside damping terms
    //
public:
    double Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed,\
        int ind_i, int ind_j, Shape_functions& Shape, double (&f)(double, double));
    double Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed,\
        Shape_functions& Shape, double(&f)(double, double));
};

class Boundary_Vector_Integral
{
    //
    // Handles integral over the boundary dS with terms in 𝜑ⱼ
    //
};

class Volume_Vector_Integral
{
    //
    // Handles integral over the volume dV with terms in 𝜑ⱼ
    // Used in the righthand-side source term
    //
public:
    double Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed,\
        int ind_i, Shape_functions& Shape, double(&f)(double, double));
    double Gaussian_Quadrature(int ind_i, std::vector<std::vector<double>>& coord_deformed,\
        Shape_functions& Shape, double(&f)(double, double));
};

class Boundary_Matrix_Integral
{
    //
    // Handles integral over the boundary dS with terms in 𝜑ᵢ×𝜑ⱼ
    //
};