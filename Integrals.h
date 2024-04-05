#pragma once
#include <cmath>
#include <iostream>
#include <vector>
#include "Shape_functions.h"
#include "Shape_fct_1D.h"
#include "Aliases.h"


class Volume_Matrix_Integral
{
    //
    // Handles integral over the volume dV with terms in a×𝜑ᵢ×𝜑ⱼ
    // Used in the left-handside damping terms
    //
public:
    double Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed,\
        int ind_i, int ind_j, Shape_functions& Shape, integrand_function f);
    virtual double Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed,\
        Shape_functions& Shape, integrand_function f);
};

class Boundary_Vector_Integral
{
    //
    // Handles integral over the boundary dS with terms in a×𝜑ⱼ
    // Used in the right-handside F terms from the natural boundary conditions
    //
public:
    double Evaluate_Integrand(double coord_master, std::vector<std::vector<double>>& coord_deformed,\
        int ind_i, Shape_fct_1D& Shape, integrand_function f);
    double Gaussian_Quadrature(int ind_i, std::vector<std::vector<double>>& coord_deformed,\
        Shape_fct_1D& Shape, integrand_function f);
};

class Volume_Vector_Integral
{
    //
    // Handles integral over the volume dV with terms in a×𝜑ⱼ
    // Used in the righthand-side source term
    //
public:
    double Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed,\
        int ind_i, Shape_functions& Shape, integrand_function f);
    double Gaussian_Quadrature(int ind_i, std::vector<std::vector<double>>& coord_deformed,\
        Shape_functions& Shape, integrand_function f);
};

class Boundary_Matrix_Integral
{
    //
    // Handles integral over the boundary dS with terms in a×𝜑ᵢ×𝜑ⱼ
    //
public:
    double Evaluate_Integrand(double coord_master, std::vector<std::vector<double>>& coord_deformed, \
        int ind_i, int ind_j, Shape_fct_1D& Shape, integrand_function f);
    double Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed, \
        Shape_fct_1D& Shape, integrand_function f);
};

class Volume_Inner_grad_Integral : public Volume_Matrix_Integral
{
    //
    // Handles integral over the volume dV with terms in a×∇𝜑ᵢ⋅∇𝜑ⱼ
    // Used in the stiffness matrix
    //
public:
    double Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed,\
        Shape_functions& Shape, integrand_function f);
};