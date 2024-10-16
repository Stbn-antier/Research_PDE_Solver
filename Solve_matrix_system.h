#pragma once
#include <iostream>
#include <vector>

class Solve_matrix_system
{
public:
	void solve_system_LU(std::vector<std::vector<double>>& A_matrix, std::vector<double>& b_vector, std::vector<double>& dP);
	void solve_system_Cholesky(std::vector<std::vector<double>>& A_matrix, std::vector<double>& b_vector, std::vector<double>& dP);
};

