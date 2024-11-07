#pragma once
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <assert.h>

class CSRMatrix;


class COOMatrix
{
public:
	std::vector<std::tuple<int, int, double>> data;
	int n_line = 0;
	int n_col = 0;
	int nnz = 0;
	bool is_sorted = 0;
	bool is_sym = false; // Only lower triangular is stored

	COOMatrix() = default;
	COOMatrix(bool is_sym);
	COOMatrix(CSRMatrix&& mat);

	void append(int i, int j, double val);
	void append(std::tuple<int, int, double>& tup) { data.push_back(tup); }
	std::tuple<int, int> shape();
	std::tuple<int, int, double> get_data(int index);
	void sum_duplicate();
	void sum_duplicate(bool col_first);
	void sort_tuple(bool col_first);
	size_t get_size() { return data.size(); }
	size_t get_nnz();
	void print_data();
	void print_data(int n_line_to_print);
	void to_dense(std::vector<std::vector<double>>& Dense);
};

class CSRMatrix
{
public:
	std::vector<int> rowptr;
	std::vector<int> colind;
	std::vector<double> values;
	int n_row = 0;
	int n_col = 0;
	int nnz = 0;

	bool is_sym = false;

	CSRMatrix() = default;
	CSRMatrix(bool is_sym);
	CSRMatrix(COOMatrix& m);
	CSRMatrix(const CSRMatrix& m);
	std::tuple<int, int> shape();
	void to_dense(std::vector<std::vector<double>>& Dense);
	void transpose(CSRMatrix& Res);
	CSRMatrix& operator= (CSRMatrix const& mat);
	void print_data();
	void print_data(int n_line);
};

class CSCMatrix {
public:
	std::vector<int> rowind;
	std::vector<int> colptr;
	std::vector<double> values;
	int n_row = 0;
	int n_col = 0;
	int nnz = 0;

	bool is_sym = false;

	CSCMatrix() = default;
	CSCMatrix(bool is_sym);
	CSCMatrix(COOMatrix& m);
	CSCMatrix(const CSCMatrix& m);
	std::tuple<int, int> shape();
	void to_dense(std::vector<std::vector<double>>& Dense);
	void transpose(CSCMatrix& Res);
	CSCMatrix& operator= (CSCMatrix const& mat);
	void print_data();
};

class DataVector {
private:
	std::vector<double> data;
public:
	int dof_per_node = 1;

	DataVector() = default;
	DataVector(int n_dof);
	DataVector(int n_nodes, int dof_node);
	DataVector(const DataVector& v);
	DataVector(std::vector<double>& v);

	void Extract_data(int dof_data, DataVector& result_vec);
	void append(double value);
	size_t size() { return data.size(); }
	double norm(); // Root sum square norm ||x||_2
	auto begin() { return data.begin(); }
	auto end() { return data.end(); }
	auto data_vec() { return data.data(); }
	auto swap(DataVector& v) { return data.swap(v.data); }
	void print();

	double operator[](int i) const { return data[i]; }
	double& operator[](int i) { return data[i]; }
	DataVector& operator= (DataVector const& v);
	DataVector& operator+= (DataVector const& v);
};

// Binary operators
// Perform operation like Mat_res = Mat1 <op> Mat2
template <class binary_op>
void CSR_binop_CSR(CSRMatrix const& mat1, CSRMatrix const& mat2, CSRMatrix& mat_res, binary_op op);

//template <typename T>
//CSRMatrix& operator+ (CSRMatrix const& mat, T const scal);
CSRMatrix operator+ (CSRMatrix& mat1, CSRMatrix& mat2);
CSRMatrix operator+ (CSRMatrix& mat1, CSRMatrix&& mat2);

//template <typename T>
//CSRMatrix& operator- (CSRMatrix const& mat, T const scal);
CSRMatrix operator- (CSRMatrix& mat1, CSRMatrix& mat2);
CSRMatrix operator- (CSRMatrix& mat1, CSRMatrix&& mat2);

template <typename T>
CSRMatrix operator* (CSRMatrix& mat, T scal);
CSRMatrix operator* (double scal, CSRMatrix& mat);

// Vector operators
DataVector operator *(CSRMatrix& mat, DataVector& vec);
DataVector operator *(CSRMatrix&& mat, DataVector& vec);
//DataVector operator *(DataVector& vec, CSRMatrix& mat);
DataVector Solve_conditioning(CSRMatrix& L, DataVector& vec);
DataVector Forward_substitution(CSRMatrix& L, DataVector& vec);
DataVector Backward_substitution(CSCMatrix& L, DataVector& vec);

template <class binary_op>
void Vec_binop_Vec(DataVector const& vec1, DataVector const& vec2, DataVector& vec_res, binary_op op);

void Vec_mult_Mat(CSRMatrix& mat, DataVector& vec, DataVector& res_vec);
void Vec_mult_Mat_sym(CSRMatrix& mat, DataVector& vec, DataVector& res_vec);

DataVector operator +(DataVector& vec1, DataVector& vec2);
DataVector operator +(DataVector& vec1, DataVector&& vec2);
DataVector operator +(DataVector&& vec1, DataVector&& vec2);
DataVector operator -(DataVector& vec1, DataVector& vec2);
DataVector operator -(DataVector& vec1, DataVector&& vec2);
DataVector operator *(DataVector& vec, double scal);
DataVector operator *(double scal, DataVector& vec);
DataVector operator *(double scal, DataVector&& vec);
double operator ^(DataVector& vec1, DataVector& vec2); // Scalar product vec1^T * vec2
double operator ^(DataVector& vec1, DataVector&& vec2); // Scalar product vec1^T * vec2

// Conversion function CSR -> CSC
void CSR_to_CSC(int n_row, int n_col, std::vector<int>& Ap, std::vector<int>& Aj, std::vector<double>& Ax\
	, std::vector<int>& Bp, std::vector<int>& Bi, std::vector<double>& Bx);
void CSR_to_CSC(CSRMatrix& A, CSCMatrix& B); // Converts CSR matrix A to CSC Matrix B, B is empty
void CSC_to_CSR(CSCMatrix& A, CSRMatrix& B);