#include "Shape_fct_3D.h"
#pragma once

void Shape_fct_3D::build_jacobian(std::vector<std::vector<double>>& coord_deformed)
{
    for (int k = 0; k < 8; k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                jacobian[k][3 * i + j] = Jacobian(gauss_points_3D[k], coord_deformed, i, j);
            }
        }
    }
}

void Shape_fct_3D::build_jacobian_det(std::vector<std::vector<double>>& coord_deformed)
{
    for (int k = 0; k < 8; k++) {
        jacobian_det[k] = \
              jacobian[k][3*0 + 0] * (jacobian[k][3*1 + 1] * jacobian[k][3*2 + 2] - jacobian[k][3*1 + 2] * jacobian[k][3*2 + 1])\
            - jacobian[k][3*0 + 1] * (jacobian[k][3*1 + 0] * jacobian[k][3*2 + 2] - jacobian[k][3*1 + 2] * jacobian[k][3*2 + 0])\
            + jacobian[k][3*0 + 2] * (jacobian[k][3*1 + 0] * jacobian[k][3*2 + 1] - jacobian[k][3*1 + 1] * jacobian[k][3*2 + 0]);
    }
}

void Shape_fct_3D::build_jacobian_invert(std::vector<std::vector<double>>& coord_deformed)
{
    for (int k = 0; k < 8; k++) {
        jacobian_invert[k][0] = +(jacobian[k][3 * 1 + 1] * jacobian[k][3 * 2 + 2] - jacobian[k][3 * 2 + 1] * jacobian[k][3 * 1 + 2]) / jacobian_det[k];
        jacobian_invert[k][1] = -(jacobian[k][3 * 1 + 0] * jacobian[k][3 * 2 + 2] - jacobian[k][3 * 2 + 0] * jacobian[k][3 * 1 + 2]) / jacobian_det[k];
        jacobian_invert[k][2] = +(jacobian[k][3 * 1 + 0] * jacobian[k][3 * 2 + 1] - jacobian[k][3 * 2 + 0] * jacobian[k][3 * 1 + 1]) / jacobian_det[k];
        jacobian_invert[k][3] = -(jacobian[k][3 * 0 + 1] * jacobian[k][3 * 2 + 2] - jacobian[k][3 * 2 + 1] * jacobian[k][3 * 0 + 2]) / jacobian_det[k];
        jacobian_invert[k][4] = +(jacobian[k][3 * 0 + 0] * jacobian[k][3 * 2 + 2] - jacobian[k][3 * 2 + 0] * jacobian[k][3 * 0 + 2]) / jacobian_det[k];
        jacobian_invert[k][5] = -(jacobian[k][3 * 0 + 0] * jacobian[k][3 * 2 + 1] - jacobian[k][3 * 2 + 0] * jacobian[k][3 * 0 + 1]) / jacobian_det[k];
        jacobian_invert[k][6] = +(jacobian[k][3 * 0 + 1] * jacobian[k][3 * 1 + 2] - jacobian[k][3 * 1 + 1] * jacobian[k][3 * 0 + 2]) / jacobian_det[k];
        jacobian_invert[k][7] = -(jacobian[k][3 * 0 + 0] * jacobian[k][3 * 1 + 2] - jacobian[k][3 * 1 + 0] * jacobian[k][3 * 0 + 2]) / jacobian_det[k];
        jacobian_invert[k][8] = +(jacobian[k][3 * 0 + 0] * jacobian[k][3 * 1 + 1] - jacobian[k][3 * 1 + 0] * jacobian[k][3 * 0 + 1]) / jacobian_det[k];
    }
}

Shape_fct_3D::Shape_fct_3D(std::vector<std::vector<double>>& coord_deformed)
{
    build_jacobian(coord_deformed);
    build_jacobian_det(coord_deformed);
    build_jacobian_invert(coord_deformed);
}

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

double Shape_fct_3D::EvaluateDerivativeDeformed(int gauss_idx, int index, int derivation_index)
{
    if (derivation_index > dim) {
        throw std::overflow_error("Index of derivation " + std::to_string(index) + " too large.");
    }
    double total = 0;
    for (int i = 0; i < dim; i++) {
        total += JacobianInvert(gauss_idx, derivation_index, i) * EvaluateDerivativeMaster(gauss_points_3D[gauss_idx], index, i);
        //total += JacobianInvert(coord_master, coord_deformed, derivation_index, i) * EvaluateDerivativeMaster(coord_master, index, i);
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

//double Shape_fct_3D::JacobianDeterminant(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed)
//{
//    return Jacobian(coord_master, coord_deformed, 0, 0) * \
//        (Jacobian(coord_master, coord_deformed, 1, 1) * Jacobian(coord_master, coord_deformed, 2, 2)\
//            - Jacobian(coord_master, coord_deformed, 1, 2) * Jacobian(coord_master, coord_deformed, 2, 1))\
//        - Jacobian(coord_master, coord_deformed, 0, 1) * \
//        (Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 2, 2)\
//            - Jacobian(coord_master, coord_deformed, 1, 2) * Jacobian(coord_master, coord_deformed, 2, 0))\
//        + Jacobian(coord_master, coord_deformed, 0, 2) * \
//        (Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 2, 1)\
//            - Jacobian(coord_master, coord_deformed, 1, 1) * Jacobian(coord_master, coord_deformed, 2, 0));
//}

//double Shape_fct_3D::JacobianInvert(std::vector<double> coord_master, std::vector<std::vector<double>> coord_deformed, int line, int column)
//{
//    if (line >= dim || column >= dim) {
//        throw std::overflow_error("Index of JacobianInvert out of range : " + std::to_string((line > column ? line : column)));
//    }
//    double determinant = JacobianDeterminant(coord_master, coord_deformed);
//    if (line == 0 && column == 0) {
//        return +(Jacobian(coord_master, coord_deformed, 1, 1) * Jacobian(coord_master, coord_deformed, 2, 2)\
//            - Jacobian(coord_master, coord_deformed, 2, 1) * Jacobian(coord_master, coord_deformed, 1, 2)) / determinant;
//    }
//    else if (line == 0 && column == 1) {
//        return -(Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 2, 2)\
//            - Jacobian(coord_master, coord_deformed, 2, 0) * Jacobian(coord_master, coord_deformed, 1, 2)) / determinant;
//    }
//    else if (line == 0 && column == 2) {
//        return +(Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 2, 1)\
//            - Jacobian(coord_master, coord_deformed, 2, 0) * Jacobian(coord_master, coord_deformed, 1, 1)) / determinant;
//    }
//    else if (line == 1 && column == 0) {
//        return -(Jacobian(coord_master, coord_deformed, 0, 1) * Jacobian(coord_master, coord_deformed, 2, 2)\
//            - Jacobian(coord_master, coord_deformed, 2, 1) * Jacobian(coord_master, coord_deformed, 0, 2)) / determinant;
//    }
//    else if (line == 1 && column == 1) {
//        return +(Jacobian(coord_master, coord_deformed, 0, 0) * Jacobian(coord_master, coord_deformed, 2, 2)\
//            - Jacobian(coord_master, coord_deformed, 2, 0) * Jacobian(coord_master, coord_deformed, 0, 2)) / determinant;
//    }
//    else if (line == 1 && column == 2) {
//        return -(Jacobian(coord_master, coord_deformed, 0, 0) * Jacobian(coord_master, coord_deformed, 2, 1)\
//            - Jacobian(coord_master, coord_deformed, 2, 0) * Jacobian(coord_master, coord_deformed, 0, 1)) / determinant;
//    }
//    else if (line == 2 && column == 0) {
//        return +(Jacobian(coord_master, coord_deformed, 0, 1) * Jacobian(coord_master, coord_deformed, 1, 2)\
//            - Jacobian(coord_master, coord_deformed, 1, 1) * Jacobian(coord_master, coord_deformed, 0, 2)) / determinant;
//    }
//    else if (line == 2 && column == 1) {
//        return -(Jacobian(coord_master, coord_deformed, 0, 0) * Jacobian(coord_master, coord_deformed, 1, 2)\
//            - Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 0, 2)) / determinant;
//    }
//    else if (line == 2 && column == 2) {
//        return +(Jacobian(coord_master, coord_deformed, 0, 0) * Jacobian(coord_master, coord_deformed, 1, 1)\
//            - Jacobian(coord_master, coord_deformed, 1, 0) * Jacobian(coord_master, coord_deformed, 0, 1)) / determinant;
//    }
//    return -1000000;
//}

double Shape_fct_3D::InnerProdGrad(std::vector<std::vector<double>>& coord_deformed, int gauss_idx, int ind_i, int ind_j, integrand_function3D f)
{
    return f(Coordinates_deformed(gauss_points_3D[gauss_idx], coord_deformed, 0), Coordinates_deformed(gauss_points_3D[gauss_idx], coord_deformed, 1), Coordinates_deformed(gauss_points_3D[gauss_idx], coord_deformed, 2))\
        * ((EvaluateDerivativeDeformed(gauss_idx, ind_i, 0) * EvaluateDerivativeDeformed(gauss_idx, ind_j, 0))\
        + (EvaluateDerivativeDeformed(gauss_idx, ind_i, 1) * EvaluateDerivativeDeformed(gauss_idx, ind_j, 1))\
        + (EvaluateDerivativeDeformed(gauss_idx, ind_i, 2) * EvaluateDerivativeDeformed(gauss_idx, ind_j, 2)))\
        * abs(JacobianDeterminant(gauss_idx));
}
