#include "Integrals.h"

const double gauss_point_dupli = 1 / sqrt(3);

std::vector<std::vector<double>> gauss_points{
{-gauss_point_dupli, -gauss_point_dupli},
{-gauss_point_dupli, gauss_point_dupli},
{gauss_point_dupli, gauss_point_dupli},
{gauss_point_dupli, -gauss_point_dupli}
};

double Volume_Matrix_Integral::Evaluate_Integrand(std::vector<double>& coord_master,\
	std::vector<std::vector<double>>& coord_deformed, int ind_i, int ind_j, Shape_functions& Shape,\
	double(&f)(double, double))
{
	return f(Shape.Coordinates_deformed(coord_master, coord_deformed, 0), Shape.Coordinates_deformed(coord_master, coord_deformed, 1))\
		* Shape.Evaluate(coord_master, ind_i) * Shape.Evaluate(coord_master, ind_j)\
		* abs(Shape.JacobianDeterminant(coord_master, coord_deformed));
}

double Volume_Matrix_Integral::Gaussian_Quadrature(int ind_i, int ind_j,\
	std::vector<std::vector<double>>& coord_deformed, Shape_functions& Shape, double(&f)(double, double))
{
	return Evaluate_Integrand(gauss_points[0], coord_deformed, ind_i, ind_j, Shape, f) +\
		Evaluate_Integrand(gauss_points[1], coord_deformed, ind_i, ind_j, Shape, f) +\
		Evaluate_Integrand(gauss_points[2], coord_deformed, ind_i, ind_j, Shape, f) +\
		Evaluate_Integrand(gauss_points[3], coord_deformed, ind_i, ind_j, Shape, f);
}

double Volume_Vector_Integral::Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed, int ind_i, Shape_functions& Shape, double(&f)(double, double))
{
	return f(Shape.Coordinates_deformed(coord_master, coord_deformed, 0), Shape.Coordinates_deformed(coord_master, coord_deformed, 1))\
		* Shape.Evaluate(coord_master, ind_i) * abs(Shape.JacobianDeterminant(coord_master, coord_deformed));
}

double Volume_Vector_Integral::Gaussian_Quadrature(int ind_i, std::vector<std::vector<double>>& coord_deformed, Shape_functions& Shape, double(&f)(double, double))
{
	return Evaluate_Integrand(gauss_points[0], coord_deformed, ind_i, Shape, f) +\
		Evaluate_Integrand(gauss_points[1], coord_deformed, ind_i, Shape, f) +\
		Evaluate_Integrand(gauss_points[2], coord_deformed, ind_i, Shape, f) +\
		Evaluate_Integrand(gauss_points[3], coord_deformed, ind_i, Shape, f);
}

double Volume_Inner_grad_Integral::Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed, Shape_functions& Shape, double(&f)(double, double))
{
	return Shape.InnerProdGrad(gauss_points[0], coord_deformed, ind_i, ind_j, f)\
		+ Shape.InnerProdGrad(gauss_points[1], coord_deformed, ind_i, ind_j, f)\
		+ Shape.InnerProdGrad(gauss_points[2], coord_deformed, ind_i, ind_j, f)\
		+ Shape.InnerProdGrad(gauss_points[3], coord_deformed, ind_i, ind_j, f);
}
