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



void Matrix_Builder::dirichlet_BC(Mesh& Reader, std::vector<std::vector<double>>& A_matrix,\
    std::vector<double>& b_vector, dirichlet_u0 u0_fct, std::vector<double>& u0_params, on_boundary on_bound)
{

    std::vector<double> T_0(Reader.num_nodes);

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

void Matrix_Builder::build_initial_T(Mesh& Reader, std::vector<double>& vect_T0, initial_T0 T0_fct)
{
    for (int k = 0; k < Reader.num_nodes; k++) {
        vect_T0[k] = T0_fct(Reader.Nodes[k]);
    }
}
