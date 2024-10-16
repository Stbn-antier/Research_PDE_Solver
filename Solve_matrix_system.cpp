#include "Solve_matrix_system.h"

void Solve_matrix_system::solve_system_LU(std::vector<std::vector<double>>& A_matrix, std::vector<double>& b_vector, std::vector<double>& dP) {
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

void Solve_matrix_system::solve_system_Cholesky(std::vector<std::vector<double>>& A_matrix, std::vector<double>& b_vector, std::vector<double>& dP)
{
    int dof = A_matrix.size();

    std::vector<double> L(dof*(dof+1)/2); // lower triangle matrix
    // Under this storage, element at index (i,j) of L is located in L[i*(i+1)/2 + j]
    std::vector<double> D(dof); // diagonal matrix
    std::vector<double> Y(dof); // auxiliary vector

    std::fill(dP.begin(), dP.end(), 0.0);

    // Initialize L:=I, diagonal equal to 1
    for (int i = 0; i < dof; i++) {
        L[i * (i + 1) / 2 + i] = 1;
    }

    // Calculate coefficients of D and L
    for (int j = 0; j < dof; j++) {
        double acc = 0;
        for (int k = 0; k < j; k++) {
            acc += L[j * (j + 1) / 2 + k] * L[j * (j + 1) / 2 + k] * D[k];
        }

        D[j] = A_matrix[j][j] - acc;

        for (int i = j + 1; i < dof; i++) {
            double acc = 0;
            for (int k = 0; k < j; k++) {
                acc += L[i * (i + 1) / 2 + k] * L[j * (j + 1) / 2 + k] * D[k];
            }

            L[i * (i + 1) / 2 + j] = (A_matrix[i][j] - acc) / D[j];
        }
    }

    // Solve first intermediate system b = L Y'
    for (int i = 0; i < dof; i++)
    {
        double acc = 0.0;
        for (int j = 0; j < i; j++) {
            acc += L[i*(i+1)/2 + j] * Y[j];
        }

        Y[i] = b_vector[i] - acc;
    }

    // Solve second intermediate system Y' = D Y
    for (int i = 0; i < dof; i++) {
        Y[i] /= D[i];
    }

    // Solve final system Y = L^T Dp
    for (int i = dof - 1; i >= 0; i--)
    {
        double acc = 0.0;
        for (int j = i; j < dof; j++) {
            acc += L[j * (j + 1) / 2 + i] * dP[j];
        }

        dP[i] = (Y[i] - acc) / L[i*(i+1)/2 + i];
    }

    L.clear();
    D.clear();
    Y.clear();
}
