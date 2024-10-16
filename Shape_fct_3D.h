#pragma once
#include <iostream>
#include <vector>
#include <string>
#include "Aliases.h"
#include "Gaussian_points.h"

class Shape_fct_3D
{
private:
    // One jacobian per gaussian point
    // Jacobian and Jacobian invert stored with indexes line, column :
    // [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2)]
    double jacobian[8][9] = { 0 };
    double jacobian_invert[8][9] = { 0 };
    double jacobian_det[8] = { 0 };

    void build_jacobian(std::vector<std::vector<double>>& coord_deformed);
    void build_jacobian_det(std::vector<std::vector<double>>& coord_deformed);
    void build_jacobian_invert(std::vector<std::vector<double>>& coord_deformed);
public:
    int dim = 3;
    int n_nodes = 8;

    Shape_fct_3D(std::vector<std::vector<double>>& coord_deformed);

    double Evaluate(std::vector<double> coord_master, int index);

    double Coordinates_deformed(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int index);

    double EvaluateDerivativeMaster(std::vector<double> coord_master, int index, int derivation_index);

    //double EvaluateDerivativeDeformed(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int index, int derivation_index);
    double EvaluateDerivativeDeformed(int gauss_idx, int index, int derivation_index);

    double Jacobian(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column);

    double Jacobian(int gauss_idx, int line, int column) { return jacobian[gauss_idx][3 * line + column]; };

    //double JacobianDeterminant(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed);

    double JacobianDeterminant(int gauss_idx) { return jacobian_det[gauss_idx]; };

    //double JacobianInvert(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column);

    double JacobianInvert(int gauss_idx, int line, int column) { return jacobian_invert[gauss_idx][3 * line + column]; };

    //double InnerProdGrad(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int ind_i, int ind_j, integrand_function3D f);
    double InnerProdGrad(std::vector<std::vector<double>>& coord_deformed, int gauss_idx, int ind_i, int ind_j, integrand_function3D f);
};