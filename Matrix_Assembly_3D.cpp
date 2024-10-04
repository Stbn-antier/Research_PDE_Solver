#include "Matrix_Assembly_3D.h"

void Matrix_Builder3D::build_matrix(Mesh& Reader, std::vector<std::vector<double>>& A, Volume_Matrix_Integral3D& Integration, integrand_function3D f)
{
    Shape_fct_3D ShapeFcts;

    for (int c = 0; c < Reader.num_Elems["hexa"]; c++) {
        // Building vector of coordinates of element
        vector<vector<double>> coords_element;
        for (int k = 0; k < ShapeFcts.n_nodes; k++) {
            coords_element.push_back(Reader.Nodes[Reader.Elems["hexa"][c].Nodes[k]]);
        }
        // Assembling the stiffness matrix for element c
        for (int i = 0; i < ShapeFcts.n_nodes; i++) {
            for (int j = 0; j < ShapeFcts.n_nodes; j++) {
                A[Reader.Elems["hexa"][c].Nodes[i]][Reader.Elems["hexa"][c].Nodes[j]] += Integration.Gaussian_Quadrature(i, j, coords_element, ShapeFcts, f);
            }
        }
        coords_element.clear();
        if ((c + 1) % (Reader.num_Elems["hexa"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["hexa"] << endl;
        }
    }
}

void Matrix_Builder3D::build_vector(Mesh& Reader, std::vector<double>& a, Volume_Vector_Integral3D& Integration, integrand_function3D f)
{
    Shape_fct_3D ShapeFcts;

    for (int c = 0; c < Reader.num_Elems["hexa"]; c++) {
        // Building vector of coordinates of element
        vector<vector<double>> coords_element;
        for (int k = 0; k < ShapeFcts.n_nodes; k++) {
            coords_element.push_back(Reader.Nodes[Reader.Elems["hexa"][c].Nodes[k]]);
        }
        // Assembling the force vector for element c
        for (int i = 0; i < ShapeFcts.n_nodes; i++) {
            a[Reader.Elems["hexa"][c].Nodes[i]] += Integration.Gaussian_Quadrature(i, coords_element, ShapeFcts, f);
        }
        coords_element.clear();
        if ((c + 1) % (Reader.num_Elems["hexa"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["hexa"] << endl;
        }
    }
}

void Matrix_Builder3D::dirichlet_BC(Mesh& Reader, std::vector<std::vector<double>>& A_matrix, std::vector<double>& b_vector, dirichlet_u0_3D u0_fct, std::vector<double>& u0_params, on_boundary on_bound)
{
    std::vector<double> T_0(Reader.num_nodes, 0.0);
    std::vector<bool> dirichlet_mask(Reader.num_nodes, 0); // Mask, if on_bound=1

    // Setting on boundaries based on expression
    for (int k = 0; k < Reader.num_nodes; k++) {
        if (on_bound(Reader, k)) {
            T_0[k] = u0_fct(Reader.Nodes[k][0], Reader.Nodes[k][1], Reader.Nodes[k][2], u0_params);
            dirichlet_mask[k] = 1;
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
        if (dirichlet_mask[k]) { // if on boundary == 1
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

void Matrix_Builder3D::neumann_BC(Mesh& Reader, std::vector<double>& f_vector, Boundary_Vector_Integral3D Integral, boundary_integrand3D f, on_boundary on_bound, std::vector<double>& params)
{
    //
    // From the mesh, assemble the neumann boundary condition as the integral on the boundary elements of dimension N-1
    //
    Shape_functions ShapeFct2D;

    for (int c = 0; c < Reader.num_Elems["quad"]; c++) {
        // If for selecting the boundary on which to apply the neumann BC
        if (on_bound(Reader, c)) {
            vector<vector<double>> coords_element;
            for (int k = 0; k < ShapeFct2D.n_nodes; k++) {
                coords_element.push_back(Reader.Nodes[Reader.Elems["quad"][c].Nodes[k]]);
            }
            for (int i = 0; i < ShapeFct2D.n_nodes; i++) {
                f_vector[Reader.Elems["quad"][c].Nodes[i]] += Integral.Gaussian_Quadrature(i, coords_element, ShapeFct2D, f, params);
            }
            coords_element.clear();
        }
        if ((c + 1) % (Reader.num_Elems["quad"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["quad"] << endl;
        }
    }
}

void Matrix_Builder3D::robin_BC(Mesh& Reader, std::vector<double>& f_vector, std::vector<std::vector<double>>& G_matrix, Boundary_Vector_Integral3D Vect_integral, Boundary_Matrix_Integral3D Mat_integral, boundary_integrand3D f_vect, boundary_integrand3D f_mat, on_boundary on_bound, std::vector<double>& params)
{
    //
    // From the mesh, assemble the Robin boundary condition as the integral on the boundary elements of dimension N-1
    // f_mat and f_vect are the integrand, respectively for the right-handside and left-handside terms
    //
    Shape_functions ShapeFct2D;

    for (int c = 0; c < Reader.num_Elems["line"]; c++) {
        // If for selecting the boundary on which to apply the neumann BC
        if (on_bound(Reader, c)) {
            vector<vector<double>> coords_element;
            for (int k = 0; k < ShapeFct2D.n_nodes; k++) {
                coords_element.push_back(Reader.Nodes[Reader.Elems["quad"][c].Nodes[k]]);
            }
            for (int i = 0; i < ShapeFct2D.n_nodes; i++) {
                f_vector[Reader.Elems["quad"][c].Nodes[i]] += Vect_integral.Gaussian_Quadrature(i, coords_element, ShapeFct2D, f_vect, params);
            }


            for (int i = 0; i < ShapeFct2D.n_nodes; i++) {
                for (int j = 0; j < ShapeFct2D.n_nodes; j++) {
                    G_matrix[Reader.Elems["quad"][c].Nodes[i]][Reader.Elems["quad"][c].Nodes[j]] += Mat_integral.Gaussian_Quadrature(i, j, coords_element, ShapeFct2D, f_mat, params);
                }
            }

            coords_element.clear();
        }
        if ((c + 1) % (Reader.num_Elems["quad"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["quad"] << endl;
        }
    }
}

void Matrix_Builder3D::build_initial_T(Mesh& Reader, std::vector<double>& vect_T0, initial_T0 T0_fct)
{
    for (int k = 0; k < Reader.num_nodes; k++) {
        vect_T0[k] = T0_fct(Reader.Nodes[k]);
    }
}
