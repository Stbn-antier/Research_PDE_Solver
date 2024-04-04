#pragma once
#include <iostream>
#include <vector>
#include "Integrals.h"
#include "Mesh.h"

class Matrix_Builder
{
private:
	int report_steps = 1; // Report progress at every 1/report_step % of progress
public:
	void build_matrix(Mesh& Reader, std::vector<std::vector<double>>& A, Volume_Matrix_Integral&, double(&f)(double, double));
	
	void build_matrix_boundary(Mesh& Reader, std::vector<std::vector<double>>& A, Boundary_Matrix_Integral& Integration, bool(&cond)(double, double));
	
	void build_vector(Mesh& Reader, std::vector<double>& a, Volume_Vector_Integral& Integration);
	
	void build_vector(Mesh& Reader, std::vector<double>& a, Boundary_Vector_Integral& Integration, bool(&cond)(double, double));
};
