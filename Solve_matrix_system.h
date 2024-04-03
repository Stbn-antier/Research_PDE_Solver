#pragma once
#include <iostream>
#include <vector>

class Solve_matrix_system
{
public:
	void solve_system(std::vector<std::vector<double>>& A_matrix, std::vector<double>& b_vector, std::vector<double>& dP);
};

