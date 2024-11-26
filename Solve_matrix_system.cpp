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

void Solve_matrix_system::sparse_Cholesky(CSRMatrix& A)
{
    CSCMatrix B;
    CSR_to_CSC(A, B);
    sparse_Cholesky(B);

    CSRMatrix temp_res;
    CSC_to_CSR(B, temp_res);
    A = temp_res;
}

void Solve_matrix_system::sparse_Cholesky(CSCMatrix& Res)
{
    assert(Res.n_col == Res.n_row); // Only works for square matrix n*n
    int n = Res.n_col; // size of square matrix

    for (int col = 0; col < n; col++) {
        int row_start = Res.colptr[col];
        int row_end = Res.colptr[col + 1];

        assert(Res.rowind[row_start] == col); // If failed, Res is not col-row sorted !

        // Divide column col by sqrt of diagonal coefficient
        double diag_coef_sqrt_invrt = 1/std::sqrt(Res.values[row_start]);
        for (int i = row_start; i < row_end; i++) {
            Res.values[i] *= diag_coef_sqrt_invrt;
        }

        // For values (i,j) where j>col, i>=j, (i,j) -= (i,col)*(j,col) if they all exist
        for (int j = col + 1; j < n; j++) {
            int j_row_start = Res.colptr[j];
            int j_row_end = Res.colptr[j + 1];
            for (int i_id = j_row_start; i_id < j_row_end; i_id++) {
                int i = Res.rowind[i_id];

                double val1 = 0;
                double val2 = 0;
                // Check through column col if (i, col) and (j, col) exist
                for (int col_i = row_start; col_i < row_end; col_i++) {
                    // Check if (i,col) exist and put value (i,col) in val1
                    if (Res.rowind[col_i] == i) {
                        val1 = Res.values[col_i];
                    }
                    // Check if (j,col) exist and put value (j,col) in val2
                    if (Res.rowind[col_i] == j) {
                        val2 = Res.values[col_i];
                    }
                }

                // Finally do (i,j) -= (i,col)*(j,col)
                Res.values[i_id] -= val1 * val2;
            }
        }
    }
}

DataVector Solve_matrix_system::PCCG(CSRMatrix& A, DataVector& b, bool log_err, int max_it, double tol)
{
    // Vector keeping error during execution
    std::vector<double> error(max_it);


    // First calculate the preconditioning L with Incomplete Cholesky
    CSRMatrix L(A);
    sparse_Cholesky(L);

    // First guess for solution x=0
    DataVector x(b.size());
    DataVector residue = b - (A * x);
    error[0] = residue.norm();
    DataVector z = Solve_conditioning(L, residue);
    DataVector p = z;

    for (int k = 0; k < max_it - 1; k++) {
        double alpha = (residue ^ z) / (p ^ (A * p));
        x += alpha * p;
        DataVector r_next = residue - alpha * (A * p);

        double err_r = r_next.norm();
        error[k + 1] = err_r;
        if (err_r < tol) {
            error.resize(k + 2);
            break;
        }

        DataVector z_next = Solve_conditioning(L, r_next);
        double beta = (r_next ^ z_next) / (residue ^ z);
        DataVector p_next = z_next + beta * p;

        residue = r_next;
        z = z_next;
        p = p_next;
    }

    if (log_err) {
        std::ofstream myfile;
        myfile.open("error.txt");
        for (int i = 0; i < error.size(); i++) {
            myfile << error[i] << std::endl;
        }
        myfile.close();
    }

    return x;
}

DataVector Solve_matrix_system::PCCG(CSRMatrix& A, DataVector& b, DataVector& x0_init, bool log_err, int max_it, double tol)
{
    // Vector keeping error during execution
    std::vector<double> error(max_it);


    // First calculate the preconditioning L with Incomplete Cholesky
    CSRMatrix L(A);
    sparse_Cholesky(L);

    // First guess for solution x=0
    DataVector x = x0_init;
    DataVector residue = b - (A * x);
    error[0] = residue.norm();
    DataVector z = Solve_conditioning(L, residue);
    DataVector p = z;

    for (int k = 0; k < max_it - 1; k++) {
        double alpha = (residue ^ z) / (p ^ (A * p));
        x += alpha * p;
        DataVector r_next = residue - alpha * (A * p);

        double err_r = r_next.norm();
        error[k + 1] = err_r;
        if (err_r < tol) {
            error.resize(k + 2);
            break;
        }

        DataVector z_next = Solve_conditioning(L, r_next);
        double beta = (r_next ^ z_next) / (residue ^ z);
        DataVector p_next = z_next + beta * p;

        residue = r_next;
        z = z_next;
        p = p_next;
    }

    if (log_err) {
        std::ofstream myfile;
        myfile.open("error.txt");
        for (int i = 0; i < error.size(); i++) {
            myfile << error[i] << std::endl;
        }
        myfile.close();
    }

    return x;
}
