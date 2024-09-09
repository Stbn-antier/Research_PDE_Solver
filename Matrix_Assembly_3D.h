#pragma once
#include <iostream>
#include <vector>
#include "Integrals_3D.h"
#include "Mesh.h"
#include "Aliases.h"

class Matrix_Builder3D
{
private:
	int report_steps = 1; // Report progress at every 1/report_step % of progress
public:
	void build_matrix(Mesh& Reader, std::vector<std::vector<double>>& A, Volume_Matrix_Integral3D& Integration, integrand_function3D f);

	void build_matrix_boundary(Mesh& Reader, std::vector<std::vector<double>>& A, Boundary_Matrix_Integral3D& Integration, on_boundary on_bound);

	void build_vector(Mesh& Reader, std::vector<double>& a, Volume_Vector_Integral3D& Integration, integrand_function3D f);

	void build_vector(Mesh& Reader, std::vector<double>& a, Boundary_Vector_Integral3D& Integration, on_boundary on_bound);

	void dirichlet_BC(Mesh& Reader, std::vector<std::vector<double>>& A_matrix, \
		std::vector<double>& b_vector, dirichlet_u0_3D u0_fct, std::vector<double>& u0_params, on_boundary on_bound);

	void neumann_BC(Mesh& Reader, std::vector<double>& f_vector, Boundary_Vector_Integral3D Integral, \
		boundary_integrand3D f, on_boundary on_bound, std::vector<double>& params);

	void robin_BC(Mesh& Reader, std::vector<double>& f_vector, std::vector<std::vector<double>>& G_matrix, \
		Boundary_Vector_Integral3D Vect_integral, Boundary_Matrix_Integral3D Mat_integral, \
		boundary_integrand3D f_vect, boundary_integrand3D f_mat, on_boundary on_bound, std::vector<double>& params);

	void build_initial_T(Mesh& Reader, std::vector<double>& vect_T0, initial_T0 T0_fct);

	void build_A_matrix();

	void build_b_vector();
};