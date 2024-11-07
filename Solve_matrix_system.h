#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include "Matrix.h"

class Solve_matrix_system
{
public:
	void solve_system_LU(std::vector<std::vector<double>>& A_matrix, std::vector<double>& b_vector, std::vector<double>& dP);
	void solve_system_Cholesky(std::vector<std::vector<double>>& A_matrix, std::vector<double>& b_vector, std::vector<double>& dP);
	void sparse_Cholesky(CSRMatrix& A_matrix);
	void sparse_Cholesky(CSCMatrix& A_matrix);
	DataVector PCCG(CSRMatrix& A, DataVector& b, bool log_err = false, int max_it = 5000, double tol = 1e-10);
};

