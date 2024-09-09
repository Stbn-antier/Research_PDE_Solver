#include "Integrals.h"

const double gauss_point = 1 / sqrt(3);

std::vector<std::vector<double>> gauss_points{
{-gauss_point, -gauss_point},
{-gauss_point, gauss_point},
{gauss_point, gauss_point},
{gauss_point, -gauss_point}
};
std::vector<double> one_d_gauss_points{
	gauss_point, -gauss_point
};

double Volume_Matrix_Integral::Evaluate_Integrand(std::vector<double>& coord_master,\
	std::vector<std::vector<double>>& coord_deformed, int ind_i, int ind_j, Shape_functions& Shape,\
	integrand_function f)
{
	return f(Shape.Coordinates_deformed(coord_master, coord_deformed, 0), Shape.Coordinates_deformed(coord_master, coord_deformed, 1))\
		* Shape.Evaluate(coord_master, ind_i) * Shape.Evaluate(coord_master, ind_j)\
		* abs(Shape.JacobianDeterminant(coord_master, coord_deformed));
}

double Volume_Matrix_Integral::Gaussian_Quadrature(int ind_i, int ind_j,\
	std::vector<std::vector<double>>& coord_deformed, Shape_functions& Shape, integrand_function f)
{
	double total = 0;
	for (int i = 0; i < Shape.n_nodes; i++) {
		total += Evaluate_Integrand(gauss_points[i], coord_deformed, ind_i, ind_j, Shape, f);
	}
	return total;
}

double Volume_Vector_Integral::Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed, int ind_i, Shape_functions& Shape, integrand_function f)
{
	return f(Shape.Coordinates_deformed(coord_master, coord_deformed, 0), Shape.Coordinates_deformed(coord_master, coord_deformed, 1))\
		* Shape.Evaluate(coord_master, ind_i) * abs(Shape.JacobianDeterminant(coord_master, coord_deformed));
}

double Volume_Vector_Integral::Gaussian_Quadrature(int ind_i, std::vector<std::vector<double>>& coord_deformed, Shape_functions& Shape, integrand_function f)
{
	//return Evaluate_Integrand(gauss_points[0], coord_deformed, ind_i, Shape, f) +\
	//	Evaluate_Integrand(gauss_points[1], coord_deformed, ind_i, Shape, f) +\
	//	Evaluate_Integrand(gauss_points[2], coord_deformed, ind_i, Shape, f) +\
	//	Evaluate_Integrand(gauss_points[3], coord_deformed, ind_i, Shape, f);
	double total = 0;
	for (int i = 0; i < Shape.n_nodes; i++) {
		total += Evaluate_Integrand(gauss_points[i], coord_deformed, ind_i, Shape, f);
	}
	return total;
}

double Volume_Inner_grad_Integral::Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed, Shape_functions& Shape, integrand_function f)
{
	//return Shape.InnerProdGrad(gauss_points[0], coord_deformed, ind_i, ind_j, f)\
	//	+ Shape.InnerProdGrad(gauss_points[1], coord_deformed, ind_i, ind_j, f)\
	//	+ Shape.InnerProdGrad(gauss_points[2], coord_deformed, ind_i, ind_j, f)\
	//	+ Shape.InnerProdGrad(gauss_points[3], coord_deformed, ind_i, ind_j, f);
	double total = 0;
	for (int i = 0; i < Shape.n_nodes; i++) {
		total += Shape.InnerProdGrad(gauss_points[i], coord_deformed, ind_i, ind_j, f);
	}
	return total;
}

double Boundary_Vector_Integral::Evaluate_Integrand(double coord_master, std::vector<std::vector<double>>& coord_deformed, int ind_i, Shape_fct_1D& Shape, boundary_integrand f, std::vector<double>& params)
{
	return f(Shape.Coordinates_deformed(coord_master, coord_deformed, 0), Shape.Coordinates_deformed(coord_master, coord_deformed, 1), params)\
		* Shape.Evaluate(coord_master, ind_i)\
		* sqrt(pow(coord_deformed[1][0] - coord_deformed[0][0], 2) + pow(coord_deformed[1][1] - coord_deformed[0][1], 2)) / 2;
}

double Boundary_Vector_Integral::Gaussian_Quadrature(int ind_i, std::vector<std::vector<double>>& coord_deformed, Shape_fct_1D& Shape, boundary_integrand f, std::vector<double>& params)
{
	return Evaluate_Integrand(one_d_gauss_points[0], coord_deformed, ind_i, Shape, f, params)\
		+ Evaluate_Integrand(one_d_gauss_points[1], coord_deformed, ind_i, Shape, f, params);
}

double Boundary_Matrix_Integral::Evaluate_Integrand(double coord_master, std::vector<std::vector<double>>& coord_deformed, int ind_i, int ind_j, Shape_fct_1D& Shape, boundary_integrand f, std::vector<double>& params)
{
	return f(Shape.Coordinates_deformed(coord_master, coord_deformed, 0), Shape.Coordinates_deformed(coord_master, coord_deformed, 1), params)\
		* Shape.Evaluate(coord_master, ind_i) * Shape.Evaluate(coord_master, ind_j)\
		* sqrt(pow(coord_deformed[1][0] - coord_deformed[0][0], 2) + pow(coord_deformed[1][1] - coord_deformed[0][1], 2)) / 2;
}

double Boundary_Matrix_Integral::Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed, Shape_fct_1D& Shape, boundary_integrand f, std::vector<double>& params)
{
	return Evaluate_Integrand(one_d_gauss_points[0], coord_deformed, ind_i, ind_j, Shape, f, params)\
		+ Evaluate_Integrand(one_d_gauss_points[1], coord_deformed, ind_i, ind_j,  Shape, f, params);
}
