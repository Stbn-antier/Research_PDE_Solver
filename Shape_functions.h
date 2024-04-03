#pragma once
#include <iostream>
#include <vector>
#include <string>

class Shape_functions
{
private:
public:
    int dim = 2;
    int n_nodes;
    int degree;

    Shape_functions(/* args */) {};
    ~Shape_functions() {};

    double Evaluate(std::vector<double> coord_master, int index);

    double EvaluateDerivativeMaster(std::vector<double> coord_master, int index, int derivation_index);

    double EvaluateDerivativeDeformed(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int index, int derivation_index);

    double Jacobian(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column);

    double JacobianDeterminant(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed);

    double JacobianInvert(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column);

    double InnerProdGrad(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int ind_i, int ind_j);

    double NormJacobianTNormal(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int normal_n);
};

