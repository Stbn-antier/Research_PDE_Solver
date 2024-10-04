#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <random>
#include <filesystem>
#include <format>
#include "Mesh.h"
#include "Shape_functions.h"
#include "Shape_fct_3D.h"
#include "Solve_matrix_system.h"
#include "Integrals_3D.h"
#include "Matrix_Assembly_3D.h"
#include "Writer.h"

using namespace std;

int n_elem = 464;
const int n_sf = 4; // Number of shape function per element
//const int n_dof = 43; // Number of nodes | à modifier à chaque fois, 505 pour t2, 143 pour t1

//const int report_step = 4; // Report progress at every 1/report_step % of progress

const double t_fin = 10; // Final time of simulation
const double n_time = 50; // Number of time steps not including t=0
const double dt = t_fin / n_time; // time step

std::vector<std::vector<double>> K; // Stiffness matrix
std::vector<double> F; // Force vector
std::vector<double> F_previous; // Force vector
std::vector<double> dP; // Solution vector
std::vector<double> dP_previous; // Solution vector at time t-1

std::vector<std::vector<double>> C; // Damping matrix for time stepping

//---------------------------------------------
// Constants for the different functions

const double k_iso = 1; // Heat transfer coefficient in W/m²
const double rho = 1; // Density of the material in kg/m^3
const double Cp = 1; // Heat capacity in J/(kg.K)

const double T_ext = 300; // Outside temperature in K
const double T_0 = 1100; // Initial temperature in K

const double a = -7;
const double b = 2;
const double c = -1;

//---------------------------------------------
const int num_loops = 100;

const double Lxmin = -0.5;
const double Lxmax = 0.5;
const double Lymin = -0.5;
const double Lymax = 0.5;
const double Lzmin = -0.5;
const double Lzmax = 0.5;

//---------------------------------------------


static double K_function(double x, double y, double z) {
    // Function in the stiffness term
    return k_iso;
}

static double C_function(double x, double y, double z) {
    return rho*Cp;
}

static double f_function(double x, double y, double z) {
    // Source due to second derivative
    // !!! Program solves −𝛥T=f, need to put minus sign
    return -2*(a+b+c);
}

static double u0_fct(double x, double y, double z, std::vector<double>& params) {
    return a*std::pow(x,2) + b*std::pow(y,2) + c*std::pow(z,2);
}

// Boundary condition functions

//static double Robin_BC_vect_fct_right(double x, double y, std::vector<double>& params) {
//    return h_coef_right(y, params) * T_ext;
//}
//
//static double Robin_BC_mat_fct_right(double x, double y, std::vector<double>& params) {
//    return h_coef_right(y, params);
//}
//
//static double Robin_BC_vect_fct_top(double x, double y, std::vector<double>& params) {
//    return h_coef_top(x, params) * T_ext;
//}
//
//static double Robin_BC_mat_fct_top(double x, double y, std::vector<double>& params) {
//    return h_coef_top(x, params);
//}
//
//static double Robin_BC_vect_fct_bottom(double x, double y, std::vector<double>& params) {
//    return hconv * T_ext;
//}
//
//static double Robin_BC_mat_fct_bottom(double x, double y, std::vector<double>& params) {
//    return hconv;
//}

// Initial time function

static double T0(std::vector<double>& coordinates) {
    return T_0;
}


//---------------------------------------------
// Functions defining where the boundaries are

//static bool left_boundary(Mesh& Reader, int index_line) {
//    // Returns TRUE if the LINE element at index index_line is on the left boundary
//    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][0] < 0 + 1e-10\
//        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][0] < 0 + 1e-10;
//}
//
//static bool right_boundary(Mesh& Reader, int index_line) {
//    // Returns TRUE if the LINE element at index index_line is on the right boundary
//    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][0] > lx2 - 1e-10\
//        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][0] > lx2 - 1e-10;
//}
//
//static bool bottom_boundary(Mesh& Reader, int index_line) {
//    // Returns TRUE if the LINE element at index index_line is on the bottom boundary
//    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][1] < 0 + 1e-10\
//        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][1] < 0 + 1e-10;
//}
//
//static bool top_boundary(Mesh& Reader, int index_line) {
//    // Returns TRUE if the LINE element at index index_line is on the top boundary
//    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][1] > Ly - 1e-10\
//        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][1] > Ly - 1e-10;
//}

static bool all_boundary(Mesh& Reader, int index_node) {
    // Returns TRUE if the NODE at index index_node in on a boundary of the domain
    return (\
        Reader.Nodes[index_node][0] < Lxmin + 1e-3 || Reader.Nodes[index_node][0] > Lxmax - 1e-3 || \
        Reader.Nodes[index_node][1] < Lymin + 1e-3 || Reader.Nodes[index_node][1] > Lymax - 1e-3 || \
        Reader.Nodes[index_node][2] < Lzmin + 1e-3 || Reader.Nodes[index_node][2] > Lzmax - 1e-3);
}

//static bool left_boundary_node(Mesh& Reader, int index_node) {
//    // Returns TRUE if the NODE at index index_node in on a boundary of the domain
//    return (Reader.Nodes[index_node][0] < 0 + 1e-10 || Reader.Nodes[index_node][1] < 0 + 1e-10);
//}

//---------------------------------------------


void build_A_matrix(std::vector<std::vector<double>>& A_matrix) {
    int len = A_matrix.size();

    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            A_matrix[i][j] = C[i][j] + dt * K[i][j] / 2;
        }
    }
}

void build_b_vector(std::vector<double>& b_vector) {
    double len = b_vector.size();

    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            b_vector[i] += (C[i][j] - dt * K[i][j] / 2) * dP_previous[j];
        }
        b_vector[i] += dt * (F_previous[i] + F[i]) / 2;
    }
}

// Functions for storing data in files
void store_1d_vector_in_file(std::string filename, vector<double>& array_to_store) {
    ofstream myfile;
    myfile.open(filename + ".txt");
    for (int i = 0; i < array_to_store.size(); i++) {
        myfile << array_to_store[i] << endl;
    }
    myfile.close();
}

void store_1d_vector_in_binary_file(std::string filename, vector<double>& array_to_store) {
    std::ofstream myfile(filename + ".bin", std::ios::binary);
    myfile.write(reinterpret_cast<const char*>(array_to_store.data()), array_to_store.size() * sizeof(double));
    
    myfile.close();
}

void store_2d_vector_in_file(std::string filename, vector<vector<double>>& array_to_store) {
    ofstream myfile;
    myfile.open(filename + ".txt");
    for (int i = 0; i < array_to_store.size(); i++) {
        for (int j = 0; j < array_to_store[0].size(); j++) {
            myfile << array_to_store[i][j] << ' ';
        }
        myfile << endl;
    }
    myfile.close();
}

int main() {

    // Read mesh in file, prepare result writer
    Mesh Reader;
    Reader.MeshReaderMSH("Mesh/mesh_patch_small.msh");

    XML_Writer Writer("results/", "Temperature");

    // Matrix incremental builder, and solver
    Matrix_Builder3D MBuild;

    Solve_matrix_system Solver;


    // Check if version is right
    cout << "File version : " << Reader.file_version << endl;
    // Check if nodes are saved well
    cout << "Number of nodes : " << Reader.num_nodes << endl;
    // for (int i=0; i<Reader.num_nodes; i++){
    //     cout << Reader.Nodes[i].get_x() << " " << Reader.Nodes[i].get_y() << " " << Reader.Nodes[i].get_z() << endl;
    // }
    // Check if elements are saved well
    cout << "Number of elements : " << Reader.num_Elems["hexa"] << endl;
    // for (int i=0; i<Reader.num_Elems["quad"]; i++){
    //     cout << Reader.Elems["quad"][i].Nodes[0] << " " << Reader.Elems["quad"][i].Nodes[1] << " " << Reader.Elems["quad"][i].Nodes[2] << " " << Reader.Elems["quad"][i].Nodes[3] << endl;
    // }
    cout << "Number of boundary elements : " << Reader.num_Elems["quad"] << endl;

    // Resize matrix to hold data
    int n_dof = Reader.num_nodes;

    K.resize(n_dof, std::vector<double>(n_dof)); // Stiffness matrix
    F.resize(n_dof); // Force vector
    F_previous.resize(n_dof); // Force vector
    dP.resize(n_dof); // Solution vector
    dP_previous.resize(n_dof); // Solution vector at time t-1
    C.resize(n_dof, std::vector<double>(n_dof)); // Damping matrix for time stepping


    // Test for Jacobian inverse and jacobian implementation
    
    Shape_fct_3D ShapeFctsTest;
    const double point = 1 / sqrt(3);
    std::vector<std::vector<double>> points_test{
    {-point, -point, -point},
    {+point, -point, -point},
    {+point, +point, -point},
    {-point, +point, -point},
    {-point, -point, +point},
    {+point, -point, +point},
    {+point, +point, +point},
    {-point, +point, +point},
    };
    const int index_point_test = 0;

    // Get nodes coordinates from element 0
    std::cout << "Tests for element 0" << std::endl;
    vector<vector<double>> coords_element;
    for (int k = 0; k < ShapeFctsTest.n_nodes; k++) {
        coords_element.push_back(Reader.Nodes[Reader.Elems["hexa"][0].Nodes[k]]);
    }
    std::cout << "Determinant" << ShapeFctsTest.JacobianDeterminant(points_test[index_point_test], coords_element) << std::endl;
    std::cout << "Jacobian values" << std::endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << "i=" << i << " j=" << j << " Value =" << ShapeFctsTest.Jacobian(points_test[index_point_test], coords_element, i, j) << std::endl;
        }
    }
    std::cout << "Inverse transpose Jacobian values" << std::endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << "i=" << i << " j=" << j << " Value =" << ShapeFctsTest.JacobianInvert(points_test[index_point_test], coords_element, i, j) << std::endl;
        }
    }


    // Building constant matrixes and vectors
    cout << "Building stiffness matrix" << endl;
    Volume_Inner_grad_Integral3D IntegrateK;
    MBuild.build_matrix(Reader, K, IntegrateK, K_function);

    cout << "Building force vector" << endl;
    Volume_Vector_Integral3D IntegrateF;
    MBuild.build_vector(Reader, F, IntegrateF, f_function);
    // Put F into previous F vector
    std::copy(F.begin(), F.end(), F_previous.begin());


    store_1d_vector_in_file("results/F", F);

    // Reset F vector before next time step
    std::fill(F.begin(), F.end(), 0.0);

    // For time dependent problems
    // Precomputing the constant matrices C and K
    cout << "Computing damping matrix" << endl;
    Volume_Matrix_Integral3D IntegrateC;
    MBuild.build_matrix(Reader, C, IntegrateC, C_function);


    store_2d_vector_in_file("results/C", C);
    store_2d_vector_in_file("results/K", K);


    // Dirichlet boundary imposition
    std::vector<double> params_dirichlet = { 0.0 };
    MBuild.dirichlet_BC(Reader, K, F_previous, u0_fct, params_dirichlet, all_boundary);


    Solver.solve_system(K, F_previous, dP);

    store_1d_vector_in_file("results/T", dP);
    store_1d_vector_in_binary_file("results/T", dP);
    Writer.Write_time_step(Reader, 1.0, "T", dP);

    //// First, we calculate the initial condition at t=0
    //MBuild.build_initial_T(Reader, dP_previous, T0);
    //store_1d_vector_in_file(path_storage + "result_step_0", dP_previous);

    // Then we loop on each timestep
    //for (int i = 0; i < n_time; i++) {
    //    double t = (i + 1) * dt;

    //    cout << "Timestep t = " << t << endl;

    //    std::vector<std::vector<double>> A_matrix(n_dof, std::vector<double>(n_dof));
    //    std::vector<double> b_vector(n_dof);
    //    std::vector<double> params_for_fct = { t , Lconv_vector[i + 1] };

    //    MBuild.build_vector(Reader, F, IntegrateF, f_function);



    //    build_A_matrix(A_matrix);
    //    build_b_vector(b_vector);

    //    // Neumann boundary imposition
    //    cout << "Neumann BC ";
    //    Boundary_Vector_Integral3D Neumann_integrator;
    //    //MBuild.neumann_BC(Reader, b_vector, Neumann_integrator, neumann_BC_fct, left_boundary, params_for_fct);
    //    //store_1d_vector_in_file("F", F);

    //    // Robin boundary imposition
    //    cout << "Robin BC ";
    //    Boundary_Matrix_Integral3D Robin_integrator;
    //    // Top boundary
    //    //MBuild.robin_BC(Reader, b_vector, A_matrix, Neumann_integrator, Robin_integrator, \
    //    //    Robin_BC_vect_fct_top, Robin_BC_mat_fct_top, top_boundary, params_for_fct);
    //    //// Right boundary
    //    //MBuild.robin_BC(Reader, b_vector, A_matrix, Neumann_integrator, Robin_integrator, \
    //    //    Robin_BC_vect_fct_right, Robin_BC_mat_fct_right, right_boundary, params_for_fct);
    //    //// Bottom boundary
    //    //MBuild.robin_BC(Reader, b_vector, A_matrix, Neumann_integrator, Robin_integrator, \
    //    //    Robin_BC_vect_fct_bottom, Robin_BC_mat_fct_bottom, bottom_boundary, params_for_fct);

    //    // Dirichlet boundary imposition
    //    //MBuild.dirichlet_BC(Reader, A_matrix, b_vector, u0_fct, params_dirichlet, left_boundary_node);

    //    params_for_fct.clear();

    //    Solver.solve_system(A_matrix, b_vector, dP);
    //    //store_1d_vector_in_file("b", b_vector);
    //    //store_2d_vector_in_file("A", A_matrix);


    //    A_matrix.clear();
    //    b_vector.clear();

    //    store_1d_vector_in_file(path_storage + "/result_step_" + to_string(i + 1), dP);
    //    Writer.Write_time_step(Reader, t, "T", dP);

    //    dP_previous.swap(dP);
    //    F_previous.swap(F);
    //    std::fill(F.begin(), F.end(), 0.0);
    //}

    return 0;
};
