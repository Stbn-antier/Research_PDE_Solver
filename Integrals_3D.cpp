#include "Integrals_3D.h"

//const double gauss_point = 1 / sqrt(3);
//
//std::vector<std::vector<double>> gauss_points_2D{
//{-gauss_point, -gauss_point},
//{-gauss_point, +gauss_point},
//{+gauss_point, +gauss_point},
//{+gauss_point, -gauss_point}
//};
//std::vector<std::vector<double>> gauss_points_3D{
//{-gauss_point, -gauss_point, -gauss_point},
//{+gauss_point, -gauss_point, -gauss_point},
//{+gauss_point, +gauss_point, -gauss_point},
//{-gauss_point, +gauss_point, -gauss_point},
//{-gauss_point, -gauss_point, +gauss_point},
//{+gauss_point, -gauss_point, +gauss_point},
//{+gauss_point, +gauss_point, +gauss_point},
//{-gauss_point, +gauss_point, +gauss_point},
//};

double Volume_Matrix_Integral3D::Evaluate_Integrand(int gauss_idx, std::vector<std::vector<double>>& coord_deformed, int ind_i, int ind_j, Shape_fct_3D& Shape, integrand_function3D f)
{
	return f(Shape.Coordinates_deformed(gauss_points_3D[gauss_idx], coord_deformed, 0), Shape.Coordinates_deformed(gauss_points_3D[gauss_idx], coord_deformed, 1), Shape.Coordinates_deformed(gauss_points_3D[gauss_idx], coord_deformed, 2))\
		* Shape.Evaluate(gauss_points_3D[gauss_idx], ind_i) * Shape.Evaluate(gauss_points_3D[gauss_idx], ind_j)\
		* abs(Shape.JacobianDeterminant(gauss_idx));
}

double Volume_Matrix_Integral3D::Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed, Shape_fct_3D& Shape, integrand_function3D f)
{
	double total = 0;
	for (int k = 0; k < Shape.n_nodes; k++) {
		total += Evaluate_Integrand(k, coord_deformed, ind_i, ind_j, Shape, f);
	}
	return total;
}

double Boundary_Vector_Integral3D::Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed, int ind_i, Shape_functions& Shape, boundary_integrand3D f, std::vector<double>& params)
{
	// Calculates the result with the area correction with ABCD the quad element in 3D, Area = ||AC x BD||/2 then divided by 4 for area of base element
	return f(Shape.Coordinates_deformed(coord_master, coord_deformed, 0), Shape.Coordinates_deformed(coord_master, coord_deformed, 1), Shape.Coordinates_deformed(coord_master, coord_deformed, 2), params)\
		* Shape.Evaluate(coord_master, ind_i)\
		* std::sqrt(\
			pow((coord_deformed[2][1] - coord_deformed[0][1]) * (coord_deformed[3][2] - coord_deformed[1][2]) - (coord_deformed[3][1] - coord_deformed[1][1]) * (coord_deformed[2][2] - coord_deformed[0][2]), 2) + \
			pow((coord_deformed[2][2] - coord_deformed[0][2]) * (coord_deformed[3][0] - coord_deformed[1][0]) - (coord_deformed[3][2] - coord_deformed[1][2]) * (coord_deformed[2][0] - coord_deformed[0][0]), 2) +\
			pow((coord_deformed[2][0] - coord_deformed[0][0]) * (coord_deformed[3][1] - coord_deformed[1][1]) - (coord_deformed[3][0] - coord_deformed[1][0]) * (coord_deformed[2][1] - coord_deformed[0][1]), 2))/8;
}

double Boundary_Vector_Integral3D::Gaussian_Quadrature(int ind_i, std::vector<std::vector<double>>& coord_deformed, Shape_functions& Shape, boundary_integrand3D f, std::vector<double>& params)
{
	double total = 0;
	for (int i = 0; i < Shape.n_nodes; i++) {
		total += Evaluate_Integrand(gauss_points_2D[i], coord_deformed, ind_i, Shape, f, params);
	}
	return total;
}

double Volume_Vector_Integral3D::Evaluate_Integrand(int gauss_idx, std::vector<std::vector<double>>& coord_deformed, int ind_i, Shape_fct_3D& Shape, integrand_function3D f)
{
	return f(Shape.Coordinates_deformed(gauss_points_3D[gauss_idx], coord_deformed, 0), Shape.Coordinates_deformed(gauss_points_3D[gauss_idx], coord_deformed, 1), Shape.Coordinates_deformed(gauss_points_3D[gauss_idx], coord_deformed, 2))\
		* Shape.Evaluate(gauss_points_3D[gauss_idx], ind_i) * abs(Shape.JacobianDeterminant(gauss_idx));
}

double Volume_Vector_Integral3D::Gaussian_Quadrature(int ind_i, std::vector<std::vector<double>>& coord_deformed, Shape_fct_3D& Shape, integrand_function3D f)
{
	double total = 0;
	for (int k = 0; k < Shape.n_nodes; k++) {
		total += Evaluate_Integrand(k, coord_deformed, ind_i, Shape, f);
	}
	return total;
}

double Boundary_Matrix_Integral3D::Evaluate_Integrand(std::vector<double>& coord_master, std::vector<std::vector<double>>& coord_deformed, int ind_i, int ind_j, Shape_functions& Shape, boundary_integrand3D f, std::vector<double>& params)
{
	return f(Shape.Coordinates_deformed(coord_master, coord_deformed, 0), Shape.Coordinates_deformed(coord_master, coord_deformed, 1), Shape.Coordinates_deformed(coord_master, coord_deformed, 2), params)\
		* Shape.Evaluate(coord_master, ind_i) * Shape.Evaluate(coord_master, ind_j)\
		* std::sqrt(\
			pow((coord_deformed[2][1] - coord_deformed[0][1]) * (coord_deformed[3][2] - coord_deformed[1][2]) - (coord_deformed[3][1] - coord_deformed[1][1]) * (coord_deformed[2][2] - coord_deformed[0][2]), 2) + \
			pow((coord_deformed[2][2] - coord_deformed[0][2]) * (coord_deformed[3][0] - coord_deformed[1][0]) - (coord_deformed[3][2] - coord_deformed[1][2]) * (coord_deformed[2][0] - coord_deformed[0][0]), 2) + \
			pow((coord_deformed[2][0] - coord_deformed[0][0]) * (coord_deformed[3][1] - coord_deformed[1][1]) - (coord_deformed[3][0] - coord_deformed[1][0]) * (coord_deformed[2][1] - coord_deformed[0][1]), 2)) / 8;
}

double Boundary_Matrix_Integral3D::Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed, Shape_functions& Shape, boundary_integrand3D f, std::vector<double>& params)
{
	double total = 0;
	for (int i = 0; i < Shape.n_nodes; i++) {
		total += Evaluate_Integrand(gauss_points_2D[i], coord_deformed, ind_i, ind_j, Shape, f, params);
	}
	return total;
}

double Volume_Inner_grad_Integral3D::Gaussian_Quadrature(int ind_i, int ind_j, std::vector<std::vector<double>>& coord_deformed, Shape_fct_3D& Shape, integrand_function3D f)
{
	double total = 0;
	for (int k = 0; k < Shape.n_nodes; k++) {
		total += Shape.InnerProdGrad(coord_deformed, k, ind_i, ind_j, f);
	}
	return total;
}
