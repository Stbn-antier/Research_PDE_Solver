#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
// #include <string>
#include "Mesh.h"
#include "Shape_functions.h"
#include "Shape_fct_1D.h"
#include "Solve_matrix_system.h"
#include "Integrals.h"

using namespace std;

const double width = 1.0;
const double height = 1.0;
const int n_width = 2; // Number of elements in x direction <-->
const int n_height = 2; // Number of elements in y direction ^|v 
int n_elem = 464;
const int n_sf = 4; // Number of shape function per element
const int n_dof = 143; // Number of nodes | à modifier à chaque fois
const double hw = 1 / width;
const double hh = 1 / height; // Finite element height&width

const int report_step = 1; // Report progress at every 1/report_step % of progress

const double t_fin = 10; // Final time of simulation
const double n_time = 10; // Number of time steps not including t=0
const double dt = t_fin / n_time; // time step

const double flux = 10.0; // Heat flux on left boundary
const double t0 = 1; // Temperature imposed on right boudary

const double gauss_point = 1 / sqrt(3); // Point to compute gaussian quadrature integration
std::vector<vector<double>> gauss_points{
    {-gauss_point, -gauss_point},
    {-gauss_point, gauss_point},
    {gauss_point, gauss_point},
    {gauss_point, -gauss_point}
};
vector<double> one_d_gauss_points{
    gauss_point, -gauss_point
};

std::vector<double> T_0(n_dof);

std::vector<std::vector<double>> K(n_dof, std::vector<double>(n_dof)); // Stiffness matrix
std::vector<double> F(n_dof); // Force vector
std::vector<double> F_previous(n_dof); // Force vector
std::vector<double> dP(n_dof); // Solution vector
std::vector<double> dP_previous(n_dof); // Solution vector at time t-1

std::vector<std::vector<double>> C(n_dof, std::vector<double>(n_dof)); // Damping matrix for time stepping

// Solving the equation KT = F with LU decomposition
// std::vector<std::vector<double>> L(n_dof, std::vector<double>(n_dof)); // lower matrix
// std::vector<std::vector<double>> U(n_dof, std::vector<double>(n_dof)); // upper matrix
// std::vector<double> Z(n_dof); // auxiliary vector

double gauss_integral_grad(int ind_i, int ind_j, std::vector<vector<double>> coord_deformed, Shape_functions Shape) {
    return Shape.InnerProdGrad(gauss_points[0], coord_deformed, ind_i, ind_j) + Shape.InnerProdGrad(gauss_points[1], coord_deformed, ind_i, ind_j)\
        + Shape.InnerProdGrad(gauss_points[2], coord_deformed, ind_i, ind_j) + Shape.InnerProdGrad(gauss_points[3], coord_deformed, ind_i, ind_j);
}

//---------------------------------------------

double C_function(double x, double y) {
    return 1.0;
}

double f_function(double x, double y) {
    return 1.0;
}

//---------------------------------------------


double alternate_function_neumann_bc(double coord_master, std::vector<vector<double>> coord_deformed, int ind_i, Shape_fct_1D Shape) {
    return 5 * Shape.Evaluate(coord_master, ind_i) * abs(coord_deformed[1][0] - coord_deformed[0][0] + coord_deformed[1][1] - coord_deformed[0][1]) / 2;
}

double alternate_gauss_integral_neumann(int ind_i, std::vector<vector<double>> coord_deformed, Shape_fct_1D Shape) {
    return alternate_function_neumann_bc(one_d_gauss_points[0], coord_deformed, ind_i, Shape)\
        + alternate_function_neumann_bc(one_d_gauss_points[1], coord_deformed, ind_i, Shape);
}

void build_K(Mesh& Reader) {
    Shape_functions ShapeFcts;
    for (int c = 0; c < Reader.num_Elems["quad"]; c++) {
        // Building vector of coordinates of element
        vector<vector<double>> coords_element;
        for (int k = 0; k < Reader.shapefct_per_node; k++) {
            coords_element.push_back(Reader.Nodes[Reader.Elems["quad"][c].Nodes[k]]);
        }
        // Assembling the stiffness matrix for element c
        for (int i = 0; i < n_sf; i++) {
            for (int j = 0; j < n_sf; j++) {
                K[Reader.Elems["quad"][c].Nodes[i]][Reader.Elems["quad"][c].Nodes[j]] += gauss_integral_grad(i, j, coords_element, ShapeFcts);
            }
        }
        coords_element.clear();
        if ((c + 1) % (Reader.num_Elems["quad"] / report_step) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["quad"] << endl;
        }
    }
}

void build_C(Mesh& Reader) {
    Volume_Matrix_Integral damping;
    Shape_functions ShapeFcts;

    for (int c = 0; c < Reader.num_Elems["quad"]; c++) {
        // Building vector of coordinates of element
        vector<vector<double>> coords_element;
        for (int k = 0; k < Reader.shapefct_per_node; k++) {
            coords_element.push_back(Reader.Nodes[Reader.Elems["quad"][c].Nodes[k]]);
        }
        // Assembling the stiffness matrix for element c
        for (int i = 0; i < n_sf; i++) {
            for (int j = 0; j < n_sf; j++) {
                C[Reader.Elems["quad"][c].Nodes[i]][Reader.Elems["quad"][c].Nodes[j]] += damping.Gaussian_Quadrature(i, j, coords_element, ShapeFcts, C_function);
            }
        }
        coords_element.clear();
        if ((c + 1) % (Reader.num_Elems["quad"] / report_step) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["quad"] << endl;
        }
    }
}

void build_F(Mesh& Reader) {
    Volume_Vector_Integral F_integral;
    Shape_functions ShapeFcts;

    for (int c = 0; c < Reader.num_Elems["quad"]; c++) {
        // Building vector of coordinates of element
        vector<vector<double>> coords_element;
        for (int k = 0; k < Reader.shapefct_per_node; k++) {
            coords_element.push_back(Reader.Nodes[Reader.Elems["quad"][c].Nodes[k]]);
        }
        // Assembling the force vector for element c
        for (int i = 0; i < n_sf; i++) {
            F[Reader.Elems["quad"][c].Nodes[i]] += F_integral.Gaussian_Quadrature(i, coords_element, ShapeFcts, f_function);
        }
        coords_element.clear();
        if ((c + 1) % (Reader.num_Elems["quad"] / report_step) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["quad"] << endl;
        }
    }
}

void neumann_BC(Mesh& Reader) {
    //
    // From the mesh, assemble the neumann boundary condition as the integral on the boundary elements of dimension N-1
    //
    Shape_fct_1D ShapeFct1D;
    for (int c = 0; c < Reader.num_Elems["line"]; c++) {
        // If for selecting the boundary on which to apply the neumann BC
        if (Reader.Nodes[Reader.Elems["line"][c].Nodes[0]][0] < 0 + 1e-10\
            && Reader.Nodes[Reader.Elems["line"][c].Nodes[1]][0] < 0 + 1e-10) {
            vector<vector<double>> coords_element;
            coords_element.push_back(Reader.Nodes[Reader.Elems["line"][c].Nodes[0]]);
            coords_element.push_back(Reader.Nodes[Reader.Elems["line"][c].Nodes[1]]);
            F[Reader.Elems["line"][c].Nodes[0]] += alternate_gauss_integral_neumann(0, coords_element, ShapeFct1D);
            F[Reader.Elems["line"][c].Nodes[1]] += alternate_gauss_integral_neumann(1, coords_element, ShapeFct1D);
            coords_element.clear();
        }
        if ((c + 1) % (Reader.num_Elems["line"] / report_step) == 0) {
            std::cout << "Progress : " << c + 1 << "/" << Reader.num_Elems["line"] << endl;
        }
    }
}

void dirichlet_BC(Mesh& Reader, double time, std::vector<std::vector<double>>& A_matrix, std::vector<double>& b_vector) {
    // Setting on the boundaries based on constant
    // for (int k=0; k < n_dof; k++){
    //     // if (Reader.Nodes[k][0] < 0 + 1e-10){
    //     //     T_0[k] = 2.0;
    //     // }
    //     if (Reader.Nodes[k][0] > 1 - 1e-10){
    //         T_0[k] = 1.0;
    //     }
    // }

    // Setting on boundaries based on expression
    for (int k = 0; k < n_dof; k++) {
        if ((Reader.Nodes[k][0] < 0 + 1e-10 || Reader.Nodes[k][0] > 1 - 1e-10 || \
            Reader.Nodes[k][1] < 0 + 1e-10 || Reader.Nodes[k][1] > 1 - 1e-10)) {
            T_0[k] = 1 + pow(Reader.Nodes[k][0], 2) / 2 + pow(Reader.Nodes[k][1], 2) + 4 * time;
        }
    }

    for (int k = 0; k < n_dof; k++) {
        // Remove K*T₀ from the F vector
        for (int i = 0; i < n_dof; i++) {
            b_vector[k] -= A_matrix[k][i] * T_0[i];
        }
    }

    for (int k = 0; k < n_dof; k++) {
        // If the dof is on the boundary, impose the temperature
        // and reset the corresponding row and column of K
        if (T_0[k] != 0.0) {
            b_vector[k] = T_0[k];

            A_matrix[k][k] = 1;
            for (int l = 0; l < n_dof; l++) {
                if (l != k) {
                    A_matrix[k][l] = 0.0;
                    A_matrix[l][k] = 0.0;
                }
            }
        }
    }
}

double Set_Initial_T(std::vector<double> coordinates) {
    return 1 + pow(coordinates[0], 2) / 2 + pow(coordinates[1], 2);
}

void build_T0(Mesh& Reader) {
    for (int k = 0; k < n_dof; k++) {
        dP_previous[k] = Set_Initial_T(Reader.Nodes[k]);
    }
}

void build_A_matrix(std::vector<std::vector<double>>& A_matrix) {
    for (int i = 0; i < n_dof; i++) {
        for (int j = 0; j < n_dof; j++) {
            A_matrix[i][j] = C[i][j] + dt * K[i][j] / 2;
        }
    }
}

void build_b_vector(std::vector<double>& b_vector) {
    for (int i = 0; i < n_dof; i++) {
        for (int j = 0; j < n_dof; j++) {
            b_vector[i] += (C[i][j] - dt * K[i][j] / 2) * dP_previous[j];
        }
        b_vector[i] += dt * (F_previous[i] + F[i]) / 2;
    }
}

// Functions for storing data in files
void store_1d_vector_in_file(std::string filename, vector<double>& array_to_store) {
    ofstream myfile;
    myfile.open("results/" + filename + ".txt");
    for (int i = 0; i < n_dof; i++) {
        myfile << array_to_store[i] << endl;
    }
    myfile.close();
}

void store_2d_vector_in_file(std::string filename, vector<vector<double>>& array_to_store) {
    ofstream myfile;
    myfile.open("results/" + filename + ".txt");
    for (int i = 0; i < n_dof; i++) {
        for (int j = 0; j < n_dof; j++) {
            myfile << array_to_store[i][j] << ' ';
        }
        myfile << endl;
    }
    myfile.close();
}

int main() {
    Mesh Reader;
    Reader.MeshReaderMSH("t1.msh");

    Solve_matrix_system Solver;

    // Check if version is right
    cout << "File version : " << Reader.file_version << endl;
    // Check if nodes are saved well
    cout << "Number of nodes : " << Reader.num_nodes << endl;
    // for (int i=0; i<Reader.num_nodes; i++){
    //     cout << Reader.Nodes[i].get_x() << " " << Reader.Nodes[i].get_y() << " " << Reader.Nodes[i].get_z() << endl;
    // }
    // Check if elements are saved well
    cout << "Number of elements : " << Reader.num_Elems["quad"] << endl;
    // for (int i=0; i<Reader.num_Elems["quad"]; i++){
    //     cout << Reader.Elems["quad"][i].Nodes[0] << " " << Reader.Elems["quad"][i].Nodes[1] << " " << Reader.Elems["quad"][i].Nodes[2] << " " << Reader.Elems["quad"][i].Nodes[3] << endl;
    // }
    cout << "Number of boundary elements : " << Reader.num_Elems["line"] << endl;
    cout << "Building stiffness matrix" << endl;

    build_K(Reader);
    // store_2d_vector_in_file("K", K);

    cout << "Building force vector" << endl;
    build_F(Reader);
    // Put F into previous F vector
    std::copy(F.begin(), F.end(), F_previous.begin());

    // cout << "Imposing Neumann BC" << endl;
    // neumann_BC(Reader);

    store_1d_vector_in_file("F", F);

    cout << "Imposing Dirichlet BC" << endl;
    // dirichlet_BC(Reader, 0, K, F);

    // Reset F vector before next time step
    std::fill(F.begin(), F.end(), 0.0);

    // For time dependent problems
    // Precomputing the constant matrices C and K
    cout << "Computing damping matrix" << endl;
    build_C(Reader);
    store_2d_vector_in_file("C", C);
    store_2d_vector_in_file("K", K);

    // First, we calculate the initial condition at t=0
    build_T0(Reader);
    store_1d_vector_in_file("result_t_0", dP_previous);

    // // Then we loop on each timestep
    for (int i = 0; i < n_time; i++) {
        double t = (i + 1) * dt;
        std::vector<std::vector<double>> A_matrix(n_dof, std::vector<double>(n_dof));
        std::vector<double> b_vector(n_dof);

        std::fill(T_0.begin(), T_0.end(), 0.0);
        build_F(Reader);

        build_A_matrix(A_matrix);
        build_b_vector(b_vector);
        dirichlet_BC(Reader, t, A_matrix, b_vector);
        Solver.solve_system(A_matrix, b_vector, dP);
        store_1d_vector_in_file("b", b_vector);
        store_2d_vector_in_file("A", A_matrix);


        A_matrix.clear();
        b_vector.clear();
        store_1d_vector_in_file("T", T_0);

        store_1d_vector_in_file("result_t_" + to_string(t), dP);
        dP_previous.swap(dP);
        F_previous.swap(F);
        std::fill(F.begin(), F.end(), 0.0);
    }

    // cout << "Solving system" << endl;
    // solve_system(K, F);

    // store_1d_vector_in_file("Results", dP);

    // store_1d_vector_in_file("T", T_0);
    return 0;
};
