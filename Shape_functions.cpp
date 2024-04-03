#include "Shape_functions.h"
#pragma once

double Shape_functions::Evaluate(std::vector<double> coord_master, int index) {
    switch (index)
    {
    case 0:
        return (1 - coord_master[0]) * (1 - coord_master[1]) / 4;
        break;
    case 1:
        return (1 + coord_master[0]) * (1 - coord_master[1]) / 4;
        break;
    case 2:
        return (1 + coord_master[0]) * (1 + coord_master[1]) / 4;
        break;
    case 3:
        return (1 - coord_master[0]) * (1 + coord_master[1]) / 4;
        break;
    default:
        throw std::overflow_error("Index of shape function " + std::to_string(index) + " too large.");
        break;
    }
}

double Shape_functions::EvaluateDerivativeMaster(std::vector<double> coord_master, int index, int derivation_index) {
    switch (derivation_index)
    {
    case 0:
        switch (index)
        {
        case 0:
            return -(1 - coord_master[1]) / 4;
            break;
        case 1:
            return (1 - coord_master[1]) / 4;
            break;
        case 2:
            return (1 + coord_master[1]) / 4;
            break;
        case 3:
            return -(1 + coord_master[1]) / 4;
            break;
        default:
            throw std::overflow_error("Index of shape function " + std::to_string(index) + " too large.");
            break;
        }
    case 1:
        switch (index)
        {
        case 0:
            return -(1 - coord_master[0]) / 4;
            break;
        case 1:
            return -(1 + coord_master[0]) / 4;
            break;
        case 2:
            return (1 + coord_master[0]) / 4;
            break;
        case 3:
            return (1 - coord_master[0]) / 4;
            break;
        default:
            throw std::overflow_error("Index of shape function " + std::to_string(index) + " too large.");
            break;
        }
    default:
        throw std::overflow_error("Index of derivation " + std::to_string(index) + " too large.");
        break;
    }
}

double Shape_functions::Jacobian(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column) {
    // Returns the i,j component of the 2x2 jacobian
    if (line >= dim || column >= dim) {
        throw std::overflow_error("Index of Jacobian out of range : " + std::to_string((line > column ? line : column)));
    }
    else if (line == 0 && column == 0) {
        double total = 0;
        for (int i = 0; i < 4; i++) {
            total += EvaluateDerivativeMaster(coord_master, i, 0) * coord_deformed[i][0];
        }
        return total;
    }
    else if (line == 1 && column == 0) {
        double total = 0;
        for (int i = 0; i < 4; i++) {
            total += EvaluateDerivativeMaster(coord_master, i, 1) * coord_deformed[i][0];
        }
        return total;
    }
    else if (line == 0 && column == 1) {
        double total = 0;
        for (int i = 0; i < 4; i++) {
            total += EvaluateDerivativeMaster(coord_master, i, 0) * coord_deformed[i][1];
        }
        return total;
    }
    else if (line == 1 && column == 1) {
        double total = 0;
        for (int i = 0; i < 4; i++) {
            total += EvaluateDerivativeMaster(coord_master, i, 1) * coord_deformed[i][1];
        }
        return total;
    }
    return -10000000;
}

double Shape_functions::JacobianDeterminant(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed) {
    return Jacobian(coord_master, coord_deformed, 0, 0) * Jacobian(coord_master, coord_deformed, 1, 1)\
        - Jacobian(coord_master, coord_deformed, 0, 1) * Jacobian(coord_master, coord_deformed, 1, 0);
}

double Shape_functions::JacobianInvert(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column) {
    // Returns the i,j component of the 2x2 jacobian
    if (line >= dim || column >= dim) {
        throw std::overflow_error("Index of JacobianInvert out of range : " + std::to_string((line > column ? line : column)));
    }
    double determinant = JacobianDeterminant(coord_master, coord_deformed);
    if (line == 0 && column == 0) {
        return Jacobian(coord_master, coord_deformed, 1, 1) / determinant;
    }
    if (line == 1 && column == 0) {
        return -Jacobian(coord_master, coord_deformed, 1, 0) / determinant;
    }
    if (line == 0 && column == 1) {
        return -Jacobian(coord_master, coord_deformed, 0, 1) / determinant;
    }
    if (line == 1 && column == 1) {
        return Jacobian(coord_master, coord_deformed, 0, 0) / determinant;
    }
    return -1000000;
}

double Shape_functions::EvaluateDerivativeDeformed(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int index, int derivation_index) {
    if (derivation_index > dim) {
        throw std::overflow_error("Index of derivation " + std::to_string(index) + " too large.");
    }
    double total = 0;
    for (int i = 0; i < dim; i++) {
        total += JacobianInvert(coord_master, coord_deformed, derivation_index, i) * EvaluateDerivativeMaster(coord_master, index, i);
    }
    return total;
}

double Shape_functions::InnerProdGrad(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int ind_i, int ind_j) {
    return (EvaluateDerivativeDeformed(coord_master, coord_deformed, ind_i, 0) * EvaluateDerivativeDeformed(coord_master, coord_deformed, ind_j, 0)\
        + EvaluateDerivativeDeformed(coord_master, coord_deformed, ind_i, 1) * EvaluateDerivativeDeformed(coord_master, coord_deformed, ind_j, 1))\
        * abs(JacobianDeterminant(coord_master, coord_deformed));
}

double Shape_functions::NormJacobianTNormal(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int normal_n) {
    switch (normal_n)
    {
    case 0:
        return abs(JacobianInvert(coord_master, coord_deformed, 1, 0) + JacobianInvert(coord_master, coord_deformed, 1, 1));
        break;
    case 1:
        return abs(JacobianInvert(coord_master, coord_deformed, 0, 0) + JacobianInvert(coord_master, coord_deformed, 0, 1));
        break;
    case 2:
        return abs(JacobianInvert(coord_master, coord_deformed, 1, 0) + JacobianInvert(coord_master, coord_deformed, 1, 1));
        break;
    case 3:
        return abs(JacobianInvert(coord_master, coord_deformed, 0, 0) + JacobianInvert(coord_master, coord_deformed, 0, 1));
        break;
    default:
        throw(std::range_error("Wrong boundary value for Normals"));
        break;
    }
    // Normal vectors for the faces of the element, with convention of numbering :
/*
    __ 2 __    ^ v
    |      |   |
    3      1   ---> u
    |      |
    __ 0 __
*/
// vector<vector<double>> normals{
//     {0, -1},
//     {1, 0},
//     {0, 1},
//     {-1, 0}};
}