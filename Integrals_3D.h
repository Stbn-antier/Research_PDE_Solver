#pragma once
#include <cmath>
#include <iostream>
#include <vector>
#include "Shape_functions.h"
#include "Shape_fct_3D.h"
#include "Aliases.h"


class Volume_Matrix_Integral3D
{
    //
    // Handles integral over the volume dV with terms in a×𝜑ᵢ×𝜑ⱼ
    // Used in the left-handside damping terms
    //
public:
    double Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed, \
        int ind_i, int ind_j, Shape_fct_3D& Shape, integrand_function3D f);
    virtual double Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed, \
        Shape_fct_3D& Shape, integrand_function3D f);
};

class Boundary_Vector_Integral3D
{
    //
    // Handles integral over the boundary dS with terms in a×𝜑ⱼ
    // Used in the right-handside F terms from the natural boundary conditions
    //
public:
    double Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed, \
        int ind_i, Shape_functions& Shape, boundary_integrand3D f, std::vector<double>& params);
    double Gaussian_Quadrature(int ind_i, std::vector<std::vector<double>>& coord_deformed, \
        Shape_functions& Shape, boundary_integrand3D f, std::vector<double>& params);
};

class Volume_Vector_Integral3D
{
    //
    // Handles integral over the volume dV with terms in a×𝜑ⱼ
    // Used in the righthand-side source term
    //
public:
    double Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed, \
        int ind_i, Shape_fct_3D& Shape, integrand_function3D f);
    double Gaussian_Quadrature(int ind_i, std::vector<std::vector<double>>& coord_deformed, \
        Shape_fct_3D& Shape, integrand_function3D f);
};

class Boundary_Matrix_Integral3D
{
    //
    // Handles integral over the boundary dS with terms in a×𝜑ᵢ×𝜑ⱼ
    //
public:
    double Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed, \
        int ind_i, int ind_j, Shape_functions& Shape, boundary_integrand3D f, std::vector<double>& params);
    double Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed, \
        Shape_functions& Shape, boundary_integrand3D f, std::vector<double>& params);
};

class Volume_Inner_grad_Integral3D : public Volume_Matrix_Integral3D
{
    //
    // Handles integral over the volume dV with terms in a×∇𝜑ᵢ⋅∇𝜑ⱼ
    // Used in the stiffness matrix
    //
public:
    double Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed, \
        Shape_fct_3D& Shape, integrand_function3D f);
};