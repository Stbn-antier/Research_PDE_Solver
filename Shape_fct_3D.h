#pragma once
#include <iostream>
#include <vector>
#include <string>
#include "Aliases.h"

class Shape_fct_3D
{
private:
public:
    int dim = 3;
    int n_nodes = 8;

    double Evaluate(std::vector<double> coord_master, int index);

    double Coordinates_deformed(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int index);

    double EvaluateDerivativeMaster(std::vector<double> coord_master, int index, int derivation_index);

    double EvaluateDerivativeDeformed(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int index, int derivation_index);

    double Jacobian(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column);

    double JacobianDeterminant(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed);

    double JacobianInvert(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column);

    double InnerProdGrad(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int ind_i, int ind_j, integrand_function3D f);
};