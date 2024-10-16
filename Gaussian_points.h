#pragma once
#include <vector>
#include <math.h>

static const double gauss_point = 1 / sqrt(3);

static std::vector<std::vector<double>> gauss_points_2D{
{-gauss_point, -gauss_point},
{-gauss_point, +gauss_point},
{+gauss_point, +gauss_point},
{+gauss_point, -gauss_point}
};
static std::vector<std::vector<double>> gauss_points_3D{
{-gauss_point, -gauss_point, -gauss_point},
{+gauss_point, -gauss_point, -gauss_point},
{+gauss_point, +gauss_point, -gauss_point},
{-gauss_point, +gauss_point, -gauss_point},
{-gauss_point, -gauss_point, +gauss_point},
{+gauss_point, -gauss_point, +gauss_point},
{+gauss_point, +gauss_point, +gauss_point},
{-gauss_point, +gauss_point, +gauss_point},
};