#pragma once
#include <iostream>
#include <vector>
#include <string>

class Shape_fct_1D
{
private:
public:
    int dim = 1;
    int n_nodes;
    int degree;

    Shape_fct_1D(/* args */) {};
    ~Shape_fct_1D() {};

    double Evaluate(double coord_master, int index);

    double EvaluateDerivativeMaster(std::vector<double> coord_master, int index, int derivation_index);

    double EvaluateDerivativeDeformed(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int index, int derivation_index);

    double Jacobian(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column);

    double JacobianDeterminant(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed);

    double JacobianInvert(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column);

    double InnerProdGrad(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int ind_i, int ind_j);

    double NormJacobianTNormal(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int normal_n);
};

