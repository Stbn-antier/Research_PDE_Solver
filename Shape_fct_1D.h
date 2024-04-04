#pragma once
#include <iostream>
#include <vector>
#include <string>

class Shape_fct_1D
{
private:
public:
    int dim = 1;
    int n_nodes = 2;
    //int degree;

    double Evaluate(double coord_master, int index);

    double Coordinates_deformed(double coord_master, std::vector<std::vector<double>>& coord_deformed, int index);
};

