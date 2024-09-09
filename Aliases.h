#pragma once
#include "Mesh.h"
// This file stores aliases for functions references in parameters


// Alias for the function a used in integrand like ∫a(...)*f(𝜑ᵢ,𝜑ⱼ,...)d𝛺
// Parameters are (x, y)
using integrand_function = double(&)(double, double);

// Alias for the function a used in integrand like ∫a(...)*f(𝜑ᵢ,𝜑ⱼ,...)d𝛺
// Parameters are (x, y, z)
using integrand_function3D = double(&)(double, double, double);

// Alias for the function a used in integrand like ∫a(...)*f(𝜑ᵢ,𝜑ⱼ,...)d𝛺
// Parameters are (x, y)
using boundary_integrand = double(&)(double, double, std::vector<double>&);

// Alias for the function a used in integrand like ∫a(...)*f(𝜑ᵢ,𝜑ⱼ,...)d𝛺
// Parameters are (x, y, z)
using boundary_integrand3D = double(&)(double, double, double, std::vector<double>&);


// Alias for the function used in Dirichlet boundary conditions u₀
// Parameters are (x, y, t)
using dirichlet_u0 = double(&)(double, double, std::vector<double>&);

// Alias for the function used in Dirichlet boundary conditions u₀
// Parameters are (x, y, t)
using dirichlet_u0_3D = double(&)(double, double, double, std::vector<double>&);


// Alias for the function used for boundary conditions to determine whether to apply the BC
using on_boundary = bool(&)(Mesh&, int);

// Alias for the functions defining the initial Temperature
using initial_T0 = double(&)(std::vector<double>&);