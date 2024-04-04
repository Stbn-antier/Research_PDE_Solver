#include "Matrix_Assembly.h"

void Matrix_Builder::build_matrix(Mesh& Reader, std::vector<std::vector<double>>& A, Volume_Matrix_Integral& Integration, double(&f)(double, double))
{
    Shape_functions ShapeFcts;

    for (int c = 0; c < Reader.num_Elems["quad"]; c++) {
        // Building vector of coordinates of element
        vector<vector<double>> coords_element;
        for (int k = 0; k < Reader.shapefct_per_node; k++) {
            coords_element.push_back(Reader.Nodes[Reader.Elems["quad"][c].Nodes[k]]);
        }
        // Assembling the stiffness matrix for element c
        for (int i = 0; i < ShapeFcts.n_nodes; i++) {
            for (int j = 0; j < ShapeFcts.n_nodes; j++) {
                A[Reader.Elems["quad"][c].Nodes[i]][Reader.Elems["quad"][c].Nodes[j]] += Integration.Gaussian_Quadrature(i, j, coords_element, ShapeFcts, f);
            }
        }
        coords_element.clear();
        if ((c + 1) % (Reader.num_Elems["quad"] / report_steps) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["quad"] << endl;
        }
    }
}

void Matrix_Builder::build_matrix_boundary(Mesh& Reader, std::vector<std::vector<double>>& A, Boundary_Matrix_Integral& Integration, bool(&cond)(double, double))
{

}
