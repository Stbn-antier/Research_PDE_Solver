#include "Matrix_Assembly.h"

void Matrix_Builder::build_matrix(Mesh& Reader, std::vector<std::vector<double>>& A, Volume_Matrix_Integral& Integration, integrand_function f)
{
    Shape_functions ShapeFcts;

    for (int c = 0; c < Reader.num_Elems["quad"]; c++) {
        // Building vector of coordinates of element
        vector<vector<double>> coords_element;
        for (int k = 0; k < Reader.shapefct_per_node; k++) {
            coords_element.push_back(Reader.Nodes[Reader.Elems["quad"][c].Nodes[k]]);
        }
        // Assembling the stiffness matrix for element c
        for (int i = 0; i < Reader.shapefct_per_node; i++) {
            for (int j = 0; j < Reader.shapefct_per_node; j++) {
                A[Reader.Elems["quad"][c].Nodes[i]][Reader.Elems["quad"][c].Nodes[j]] += Integration.Gaussian_Quadrature(i, j, coords_element, ShapeFcts, f);
            }
        }
        coords_element.clear();
        if ((c + 1) % (Reader.num_Elems["quad"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["quad"] << endl;
        }
    }
}

void Matrix_Builder::build_matrix(Mesh& Reader, COOMatrix& A, Volume_Matrix_Integral& Integration, integrand_function f, double tol)
{
    assert(A.is_sym); // Check if A is built as a symmetric matrix

    for (int c = 0; c < Reader.num_Elems["quad"]; c++) {
        // Building vector of coordinates of element
        vector<vector<double>> coords_element;
        for (int k = 0; k < Reader.node_per_quad; k++) {
            coords_element.push_back(Reader.Nodes[Reader.Elems["quad"][c].Nodes[k]]);
        }

        Shape_functions ShapeFcts;

        // Assembling the stiffness matrix for element c
        for (int i = 0; i < ShapeFcts.n_nodes; i++) {
            for (int j = 0; j < ShapeFcts.n_nodes; j++) {
                int index_i = Reader.Elems["quad"][c].Nodes[i];
                int index_j = Reader.Elems["quad"][c].Nodes[j];
                if (index_i >= index_j) {
                    double value = Integration.Gaussian_Quadrature(i, j, coords_element, ShapeFcts, f);
                    if (std::abs(value) > tol) {
                        A.append(index_i, index_j, value);
                    }
                }
            }
        }
        coords_element.clear();
        /*if ((c + 1) % (Reader.num_Elems["hexa"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["hexa"] << endl;
        }*/
    }
}

void Matrix_Builder::build_matrix_boundary(Mesh& Reader, std::vector<std::vector<double>>& A, Boundary_Matrix_Integral& Integration, on_boundary on_bound)
{

}

void Matrix_Builder::build_vector(Mesh& Reader, std::vector<double>& a, Volume_Vector_Integral& Integration, integrand_function f)
{
    Shape_functions ShapeFcts;

    for (int c = 0; c < Reader.num_Elems["quad"]; c++) {
        // Building vector of coordinates of element
        vector<vector<double>> coords_element;
        for (int k = 0; k < Reader.shapefct_per_node; k++) {
            coords_element.push_back(Reader.Nodes[Reader.Elems["quad"][c].Nodes[k]]);
        }
        // Assembling the force vector for element c
        for (int i = 0; i < Reader.shapefct_per_node; i++) {
            a[Reader.Elems["quad"][c].Nodes[i]] += Integration.Gaussian_Quadrature(i, coords_element, ShapeFcts, f);
        }
        coords_element.clear();
        if ((c + 1) % (Reader.num_Elems["quad"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["quad"] << endl;
        }
    }
}

void Matrix_Builder::build_vector(Mesh& Reader, DataVector& a, Volume_Vector_Integral& Integration, integrand_function f)
{
    for (int c = 0; c < Reader.num_Elems["quad"]; c++) {
        // Building vector of coordinates of element
        vector<vector<double>> coords_element;
        for (int k = 0; k < Reader.node_per_quad; k++) {
            coords_element.push_back(Reader.Nodes[Reader.Elems["quad"][c].Nodes[k]]);
        }

        Shape_functions ShapeFcts;

        // Assembling the force vector for element c
        for (int i = 0; i < ShapeFcts.n_nodes; i++) {
            a[Reader.Elems["quad"][c].Nodes[i]] += Integration.Gaussian_Quadrature(i, coords_element, ShapeFcts, f);
        }
        coords_element.clear();
        /*if ((c + 1) % (Reader.num_Elems["hexa"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["hexa"] << endl;
        }*/
    }
}



void Matrix_Builder::dirichlet_BC(Mesh& Reader, std::vector<std::vector<double>>& A_matrix,\
    std::vector<double>& b_vector, dirichlet_u0 u0_fct, std::vector<double>& u0_params, on_boundary on_bound)
{

    std::vector<double> T_0(Reader.num_nodes, 0.0);

    // Setting on boundaries based on expression
    for (int k = 0; k < Reader.num_nodes; k++) {
        if (on_bound(Reader, k)) {
            T_0[k] = u0_fct(Reader.Nodes[k][0], Reader.Nodes[k][1], u0_params);
        }
    }

    for (int k = 0; k < Reader.num_nodes; k++) {
        // Remove K*T₀ from the b vector
        for (int i = 0; i < Reader.num_nodes; i++) {
            b_vector[k] -= A_matrix[k][i] * T_0[i];
        }
    }

    for (int k = 0; k < Reader.num_nodes; k++) {
        // If the dof is on the boundary, impose the temperature
        // and reset the corresponding row and column of K
        if (T_0[k] != 0.0) {
            b_vector[k] = T_0[k];

            A_matrix[k][k] = 1;
            for (int l = 0; l < Reader.num_nodes; l++) {
                if (l != k) {
                    A_matrix[k][l] = 0.0;
                    A_matrix[l][k] = 0.0;
                }
            }
        }
    }

    T_0.clear();
}

void Matrix_Builder::neumann_BC(Mesh& Reader, std::vector<double>& f_vector, Boundary_Vector_Integral Integral, boundary_integrand f, on_boundary on_bound, std::vector<double>& params)
{
    //
    // From the mesh, assemble the neumann boundary condition as the integral on the boundary elements of dimension N-1
    //
    Shape_fct_1D ShapeFct1D;

    for (int c = 0; c < Reader.num_Elems["line"]; c++) {
        // If for selecting the boundary on which to apply the neumann BC
        if (on_bound(Reader, c)) {
            vector<vector<double>> coords_element;
            coords_element.push_back(Reader.Nodes[Reader.Elems["line"][c].Nodes[0]]);
            coords_element.push_back(Reader.Nodes[Reader.Elems["line"][c].Nodes[1]]);
            f_vector[Reader.Elems["line"][c].Nodes[0]] += Integral.Gaussian_Quadrature(0, coords_element, ShapeFct1D, f, params);
            f_vector[Reader.Elems["line"][c].Nodes[1]] += Integral.Gaussian_Quadrature(1, coords_element, ShapeFct1D, f, params);
            coords_element.clear();
        }
        if ((c + 1) % (Reader.num_Elems["line"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["line"] << endl;
        }
    }
}

void Matrix_Builder::robin_BC(Mesh& Reader, std::vector<double>& f_vector, std::vector<std::vector<double>>& G_matrix, \
    Boundary_Vector_Integral Vect_integral, Boundary_Matrix_Integral Mat_integral, \
    boundary_integrand f_vect, boundary_integrand f_mat, on_boundary on_bound, std::vector<double>& params)
{
    //
    // From the mesh, assemble the Robin boundary condition as the integral on the boundary elements of dimension N-1
    // f_mat and f_vect are the integrand, respectively for the right-handside and left-handside terms
    //
    Shape_fct_1D ShapeFct1D;

    for (int c = 0; c < Reader.num_Elems["line"]; c++) {
        // If for selecting the boundary on which to apply the neumann BC
        if (on_bound(Reader, c)) {
            vector<vector<double>> coords_element;
            coords_element.push_back(Reader.Nodes[Reader.Elems["line"][c].Nodes[0]]);
            coords_element.push_back(Reader.Nodes[Reader.Elems["line"][c].Nodes[1]]);
            f_vector[Reader.Elems["line"][c].Nodes[0]] += Vect_integral.Gaussian_Quadrature(0, coords_element, ShapeFct1D, f_vect, params);
            f_vector[Reader.Elems["line"][c].Nodes[1]] += Vect_integral.Gaussian_Quadrature(1, coords_element, ShapeFct1D, f_vect, params);
            

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    G_matrix[Reader.Elems["line"][c].Nodes[i]][Reader.Elems["line"][c].Nodes[j]] += Mat_integral.Gaussian_Quadrature(i, j, coords_element, ShapeFct1D, f_mat, params);
                }
            }

            coords_element.clear();
        }
        if ((c + 1) % (Reader.num_Elems["line"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["line"] << endl;
        }
    }
}

void Matrix_Builder::robin_BC(Mesh& Reader, DataVector& f_vector, COOMatrix& G_matrix, Boundary_Vector_Integral& Vect_integral, Boundary_Matrix_Integral& Mat_integral, boundary_integrand f_vect, boundary_integrand f_mat, on_boundary on_bound, std::vector<double>& params, double tol)
{
    assert(G_matrix.is_sym); // Check if A is built as a symmetric matrix
    if (not G_matrix.is_sym) {
        std::cout << "Issue G not symmetric" << std::endl;
    }
    //
    // From the mesh, assemble the Robin boundary condition as the integral on the boundary elements of dimension N-1
    // f_mat and f_vect are the integrand, respectively for the right-handside and left-handside terms
    //
    Shape_fct_1D ShapeFct1D;

    for (int c = 0; c < Reader.num_Elems["line"]; c++) {
        // If for selecting the boundary on which to apply the neumann BC
        if (on_bound(Reader, c)) {
            vector<vector<double>> coords_element;
            for (int k = 0; k < ShapeFct1D.n_nodes; k++) {
                coords_element.push_back(Reader.Nodes[Reader.Elems["line"][c].Nodes[k]]);
            }
            for (int i = 0; i < ShapeFct1D.n_nodes; i++) {
                f_vector[Reader.Elems["line"][c].Nodes[i]] += Vect_integral.Gaussian_Quadrature(i, coords_element, ShapeFct1D, f_vect, params);
            }


            for (int i = 0; i < ShapeFct1D.n_nodes; i++) {
                for (int j = 0; j < ShapeFct1D.n_nodes; j++) {
                    int index_i = Reader.Elems["line"][c].Nodes[i];
                    int index_j = Reader.Elems["line"][c].Nodes[j];
                    if (index_i >= index_j) {
                        double value = Mat_integral.Gaussian_Quadrature(i, j, coords_element, ShapeFct1D, f_mat, params);
                        if (std::abs(value) > tol) {
                            G_matrix.append(index_i, index_j, value);
                        }
                    }
                }
            }

            coords_element.clear();
        }
        /*if ((c + 1) % (Reader.num_Elems["quad"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["quad"] << endl;
        }*/
    }
}

void Matrix_Builder::build_initial_T(Mesh& Reader, std::vector<double>& vect_T0, initial_T0 T0_fct)
{
    for (int k = 0; k < Reader.num_nodes; k++) {
        vect_T0[k] = T0_fct(Reader.Nodes[k]);
    }
}

void Matrix_Builder::build_initial_T(Mesh& Reader, DataVector& vect_T0, initial_T0 T0_fct)
{
    for (int k = 0; k < Reader.num_nodes; k++) {
        vect_T0[k] = T0_fct(Reader.Nodes[k]);
    }
}
