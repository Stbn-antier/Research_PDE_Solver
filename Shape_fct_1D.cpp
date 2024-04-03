#include "Shape_fct_1D.h"

double Shape_fct_1D::Evaluate(double coord_master, int index) {
    switch (index)
    {
    case 0:
        return (1 - coord_master) / 2;
        break;
    case 1:
        return (1 + coord_master) / 2;
        break;
    default:
        throw std::overflow_error("Index of shape function " + std::to_string(index) + " too large.");
        break;
    }
}

double Shape_fct_1D::Coordinates_deformed(double coord_master, std::vector<std::vector<double>>& coord_deformed, int index)
{
    // Returns the position of the x_i coordinate in the deformed element from the nodes and reference element coordinate ksi
    double position = 0.0;
    for (int k = 0; k < n_nodes; k++) {
        position += coord_deformed[k][index] * Evaluate(coord_master, k);
    }
    return position;
}
