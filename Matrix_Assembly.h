#pragma once
#include <iostream>
#include <vector>
#include "Integrals.h"
#include "Mesh.h"
#include "Aliases.h"

class Matrix_Builder
{
private:
	int report_steps = 1; // Report progress at every 1/report_step % of progress
public:
	void build_matrix(Mesh& Reader, std::vector<std::vector<double>>& A, Volume_Matrix_Integral& Integration, integrand_function f);
	
	void build_matrix_boundary(Mesh& Reader, std::vector<std::vector<double>>& A, Boundary_Matrix_Integral& Integration, on_boundary on_bound);
	
	void build_vector(Mesh& Reader, std::vector<double>& a, Volume_Vector_Integral& Integration, integrand_function f);
	
	void build_vector(Mesh& Reader, std::vector<double>& a, Boundary_Vector_Integral& Integration, on_boundary on_bound);

	void dirichlet_BC(Mesh& Reader, std::vector<std::vector<double>>& A_matrix, \
		std::vector<double>& b_vector, dirichlet_u0 u0_fct, std::vector<double>& u0_params, on_boundary on_bound);

	void build_initial_T(Mesh& Reader, std::vector<double>& vect_T0, initial_T0 T0_fct);

	void build_A_matrix();

	void build_b_vector();
};
