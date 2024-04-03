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