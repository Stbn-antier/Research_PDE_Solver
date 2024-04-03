#include "Solve_matrix_system.h"

void Solve_matrix_system::solve_system(std::vector<std::vector<double>>& A_matrix, std::vector<double>& b_vector, std::vector<double>& dP) {
    int dof = A_matrix.size();

    std::vector<std::vector<double>> L(dof, std::vector<double>(dof)); // lower matrix
    std::vector<std::vector<double>> U(dof, std::vector<double>(dof)); // upper matrix
    std::vector<double> Z(dof); // auxiliary vector

    std::fill(dP.begin(), dP.end(), 0.0);

    // init L := I
    for (int i = 0; i < dof; i++) {
        L[i][i] = 1.0;
    }

    // find L, U where L * U = A_matrix
    for (int i1 = 0; i1 < dof; i1++) {
        double acc = 0.0;
        for (int i2 = 0; i2 < dof; i2++)
            acc += L[i1][i2] * U[i2][i1];
        U[i1][i1] = A_matrix[i1][i1] - acc;

        for (int i2 = i1 + 1; i2 < dof; i2++) {
            acc = 0.0;
            for (int i3 = 0; i3 < i1; i3++)
                acc += L[i1][i3] * U[i3][i2];
            U[i1][i2] = A_matrix[i1][i2] - acc;

            acc = 0.0;
            for (int i3 = 0; i3 < i1; i3++)
                acc += L[i2][i3] * U[i3][i1];
            L[i2][i1] = (A_matrix[i2][i1] - acc) / U[i1][i1];
        }
    }

    // finally find result
    for (int i1 = 0; i1 < dof; i1++)
    {
        // find Z where L * Z = b_vector
        double acc = 0.0;
        for (int i2 = 0; i2 < i1; i2++)
            acc += L[i1][i2] * Z[i2];
        Z[i1] = b_vector[i1] - acc;
    }
    for (int i1 = dof - 1; i1 >= 0; i1--)
    {
        // find dP where U * dP = Z
        double acc = 0.0;
        for (int i2 = i1; i2 < dof; i2++)
            acc += U[i1][i2] * dP[i2];
        dP[i1] = (Z[i1] - acc) / U[i1][i1];
    }

    L.clear();
    U.clear();
    Z.clear();
}