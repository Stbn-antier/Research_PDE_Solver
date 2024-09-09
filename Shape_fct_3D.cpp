#include "Shape_fct_3D.h"
#pragma once

double Shape_fct_3D::Evaluate(std::vector<double> coord_master, int index)
{
    switch (index)
    {
    case 0:
        return (1 - coord_master[0]) * (1 - coord_master[1]) * (1 - coord_master[2]) / 8;
        break;
    case 1:
        return (1 + coord_master[0]) * (1 - coord_master[1]) * (1 - coord_master[2]) / 8;
        break;
    case 2:
        return (1 + coord_master[0]) * (1 + coord_master[1]) * (1 - coord_master[2]) / 8;
        break;
    case 3:
        return (1 - coord_master[0]) * (1 + coord_master[1]) * (1 - coord_master[2]) / 8;
        break;
    case 4:
        return (1 - coord_master[0]) * (1 - coord_master[1]) * (1 + coord_master[2]) / 8;
        break;
    case 5:
        return (1 + coord_master[0]) * (1 - coord_master[1]) * (1 + coord_master[2]) / 8;
        break;
    case 6:
        return (1 + coord_master[0]) * (1 + coord_master[1]) * (1 + coord_master[2]) / 8;
        break;
    case 7:
        return (1 - coord_master[0]) * (1 + coord_master[1]) * (1 + coord_master[2]) / 8;
        break;
    default:
        throw std::overflow_error("Index of shape function " + std::to_string(index) + " too large.");
        break;
    }
}

double Shape_fct_3D::Coordinates_deformed(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int index)
{
    // Returns the position of the x_i coordinate in the deformed element from the nodes and reference element coordinates
    double position = 0.0;
    for (int k = 0; k < n_nodes; k++) {
        position += coord_deformed[k][index] * Evaluate(coord_master, k);
    }
    return position;
}

double Shape_fct_3D::EvaluateDerivativeMaster(std::vector<double> coord_master, int index, int derivation_index)
{
    switch (derivation_index)
    {
    case 0:
        switch (index)
        {
        case 0:
            return -(1 - coord_master[1]) * (1 - coord_master[2]) / 8;
            break;
        case 1:
            return (1 - coord_master[1]) * (1 - coord_master[2]) / 8;
            break;
        case 2:
            return (1 + coord_master[1]) * (1 - coord_master[2]) / 8;
            break;
        case 3:
            return -(1 + coord_master[1]) * (1 - coord_master[2]) / 8;
            break;
        case 4:
            return -(1 - coord_master[1]) * (1 + coord_master[2]) / 8;
            break;
        case 5:
            return (1 - coord_master[1]) * (1 + coord_master[2]) / 8;
            break;
        case 6:
            return (1 + coord_master[1]) * (1 + coord_master[2]) / 8;
            break;
        case 7:
            return -(1 + coord_master[1]) * (1 + coord_master[2]) / 8;
            break;
        default:
            throw std::overflow_error("Index of shape function " + std::to_string(index) + " too large.");
            break;
        }
    case 1:
        switch (index)
        {
        case 0:
            return -(1 - coord_master[0]) * (1 - coord_master[2]) / 8;
            break;
        case 1:
            return -(1 + coord_master[0]) * (1 - coord_master[2]) / 8;
            break;
        case 2:
            return (1 + coord_master[0]) * (1 - coord_master[2]) / 8;
            break;
        case 3:
            return (1 - coord_master[0]) * (1 - coord_master[2]) / 8;
            break;
        case 4:
            return -(1 - coord_master[0]) * (1 + coord_master[2]) / 8;
            break;
        case 5:
            return -(1 + coord_master[0]) * (1 + coord_master[2]) / 8;
            break;
        case 6:
            return (1 + coord_master[0]) * (1 + coord_master[2]) / 8;
            break;
        case 7:
            return (1 - coord_master[0]) * (1 + coord_master[2]) / 8;
            break;
        default:
            throw std::overflow_error("Index of shape function " + std::to_string(index) + " too large.");
            break;
        }
    case 2:
        switch (index)
        {
        case 0:
            return -(1 - coord_master[0]) * (1 - coord_master[1]) / 8;
            break;
        case 1:
            return -(1 + coord_master[0]) * (1 - coord_master[1]) / 8;
            break;
        case 2:
            return -(1 + coord_master[0]) * (1 + coord_master[1]) / 8;
            break;
        case 3:
            return -(1 - coord_master[0]) * (1 + coord_master[1]) / 8;
            break;
        case 4:
            return (1 - coord_master[0]) * (1 - coord_master[1]) / 8;
            break;
        case 5:
            return (1 + coord_master[0]) * (1 - coord_master[1]) / 8;
            break;
        case 6:
            return (1 + coord_master[0]) * (1 + coord_master[1]) / 8;
            break;
        case 7:
            return (1 - coord_master[0]) * (1 + coord_master[1]) / 8;
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

double Shape_fct_3D::EvaluateDerivativeDeformed(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int index, int derivation_index)
{
    if (derivation_index > dim) {
        throw std::overflow_error("Index of derivation " + std::to_string(index) + " too large.");
    }
    double total = 0;
    for (int i = 0; i < dim; i++) {
        total += JacobianInvert(coord_master, coord_deformed, derivation_index, i) * EvaluateDerivativeMaster(coord_master, index, i);
    }
    return total;
}

double Shape_fct_3D::Jacobian(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column)
{
    // Returns the i,j component of the 3x3 jacobian
    double total = 0;
    for (int i = 0; i < n_nodes; i++) {
        total += EvaluateDerivativeMaster(coord_master, i, line) * coord_deformed[i][column];
    }
    return total;
}

double Shape_fct_3D::JacobianDeterminant(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed)
{
    return Jacobian(coord_master, coord_deformed, 0, 0) * \
        (Jacobian(coord_master, coord_deformed, 1, 1) * Jacobian(coord_master, coord_deformed, 2, 2)\
            - Jacobian(coord_master, coord_deformed, 1, 2) * Jacobian(coord_master, coord_deformed, 2, 1))\
        - Jacobian(coord_master, coord_deformed, 0, 1) * \
        (Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 2, 2)\
            - Jacobian(coord_master, coord_deformed, 1, 2) * Jacobian(coord_master, coord_deformed, 2, 0))\
        + Jacobian(coord_master, coord_deformed, 0, 2) * \
        (Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 2, 1)\
            - Jacobian(coord_master, coord_deformed, 1, 1) * Jacobian(coord_master, coord_deformed, 2, 0));
}

double Shape_fct_3D::JacobianInvert(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column)
{
    if (line >= dim || column >= dim) {
        throw std::overflow_error("Index of JacobianInvert out of range : " + std::to_string((line > column ? line : column)));
    }
    double determinant = JacobianDeterminant(coord_master, coord_deformed);
    if (line == 0 && column == 0) {
        return +(Jacobian(coord_master, coord_deformed, 1, 1) * Jacobian(coord_master, coord_deformed, 2, 2)\
            - Jacobian(coord_master, coord_deformed, 2, 1) * Jacobian(coord_master, coord_deformed, 1, 2)) / determinant;
    }
    else if (line == 0 && column == 1) {
        return -(Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 2, 2)\
            - Jacobian(coord_master, coord_deformed, 2, 0) * Jacobian(coord_master, coord_deformed, 1, 2)) / determinant;
    }
    else if (line == 0 && column == 2) {
        return +(Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 2, 1)\
            - Jacobian(coord_master, coord_deformed, 2, 0) * Jacobian(coord_master, coord_deformed, 1, 1)) / determinant;
    }
    else if (line == 1 && column == 0) {
        return -(Jacobian(coord_master, coord_deformed, 0, 1) * Jacobian(coord_master, coord_deformed, 2, 2)\
            - Jacobian(coord_master, coord_deformed, 2, 1) * Jacobian(coord_master, coord_deformed, 0, 2)) / determinant;
    }
    else if (line == 1 && column == 1) {
        return +(Jacobian(coord_master, coord_deformed, 0, 0) * Jacobian(coord_master, coord_deformed, 2, 2)\
            - Jacobian(coord_master, coord_deformed, 2, 0) * Jacobian(coord_master, coord_deformed, 0, 2)) / determinant;
    }
    else if (line == 1 && column == 2) {
        return -(Jacobian(coord_master, coord_deformed, 0, 0) * Jacobian(coord_master, coord_deformed, 2, 1)\
            - Jacobian(coord_master, coord_deformed, 2, 0) * Jacobian(coord_master, coord_deformed, 0, 1)) / determinant;
    }
    else if (line == 2 && column == 0) {
        return +(Jacobian(coord_master, coord_deformed, 0, 1) * Jacobian(coord_master, coord_deformed, 1, 2)\
            - Jacobian(coord_master, coord_deformed, 1, 1) * Jacobian(coord_master, coord_deformed, 0, 2)) / determinant;
    }
    else if (line == 2 && column == 1) {
        return -(Jacobian(coord_master, coord_deformed, 0, 0) * Jacobian(coord_master, coord_deformed, 1, 2)\
            - Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 0, 2)) / determinant;
    }
    else if (line == 2 && column == 2) {
        return +(Jacobian(coord_master, coord_deformed, 0, 0) * Jacobian(coord_master, coord_deformed, 1, 1)\
            - Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 0, 1)) / determinant;
    }
    return -1000000;
}

double Shape_fct_3D::InnerProdGrad(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int ind_i, int ind_j, integrand_function3D f)
{
    return f(Coordinates_deformed(coord_master, coord_deformed, 0), Coordinates_deformed(coord_master, coord_deformed, 1), Coordinates_deformed(coord_master, coord_deformed, 2))\
        * (EvaluateDerivativeDeformed(coord_master, coord_deformed, ind_i, 0) * EvaluateDerivativeDeformed(coord_master, coord_deformed, ind_j, 0)\
        + EvaluateDerivativeDeformed(coord_master, coord_deformed, ind_i, 1) * EvaluateDerivativeDeformed(coord_master, coord_deformed, ind_j, 1)\
        + EvaluateDerivativeDeformed(coord_master, coord_deformed, ind_i, 2) * EvaluateDerivativeDeformed(coord_master, coord_deformed, ind_j, 2))\
        * abs(JacobianDeterminant(coord_master, coord_deformed));
}
