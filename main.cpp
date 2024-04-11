#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include "Mesh.h"
#include "Shape_functions.h"
#include "Shape_fct_1D.h"
#include "Solve_matrix_system.h"
#include "Integrals.h"
#include "Matrix_Assembly.h"

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

std::vector<std::vector<double>> K(n_dof, std::vector<double>(n_dof)); // Stiffness matrix
std::vector<double> F(n_dof); // Force vector
std::vector<double> F_previous(n_dof); // Force vector
std::vector<double> dP(n_dof); // Solution vector
std::vector<double> dP_previous(n_dof); // Solution vector at time t-1

std::vector<std::vector<double>> C(n_dof, std::vector<double>(n_dof)); // Damping matrix for time stepping

//---------------------------------------------
// Constants for the different functions

const double k_iso = 50; // Heat transfer coefficient in W/m²
const double rho = 1; // Density of the material in kg/m^3
const double Cp = 500; // Heat capacity in J/(kg.K)
const double q_coef = 5000; // Inward heat flux in W/m²
const double T_ext = 200; // Outside temperature in K
const double T_0 = 300; // Initial temperature in K

const double h_coef = 10; // Heat transfer coefficient for the right boundary in W/(m².K)

//---------------------------------------------


static double K_function(double x, double y) {
    // Function in the stiffness term
    return k_iso;
}

static double C_function(double x, double y) {
    return rho*Cp;
}

static double f_function(double x, double y) {
    // No source terms
    return 0.0;
}

static double neumann_BC_fct(double x, double y, std::vector<double>& params) {
    // Function in the integrand for the Neumann BC
    if (y < 0.7 && y > 0.3 && params[0] >= 3) {
        return q_coef;
    };
    return 0;
}

static double Robin_BC_vect_fct(double x, double y, std::vector<double>& params) {
    // Function in the integrand for the Neumann BC
    return h_coef * T_ext;
}

static double Robin_BC_mat_fct(double x, double y, std::vector<double>& params) {
    // No source terms
    return h_coef;
}

//static double u0_fct(double x, double y, std::vector<double>& params) {
//    // Function u0 in the Dirichlet boundary condition : T=u₀ on 𝛤
//    return 1 + pow(x, 2) / 2 + pow(y, 2) + 4 * params[0];
//}

static double T0(std::vector<double>& coordinates) {
    return T_0;
}


//---------------------------------------------

static bool left_boundary(Mesh& Reader, int index_line) {
    // Returns TRUE if the LINE element at index index_line is on the left boundary
    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][0] < 0 + 1e-10\
        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][0] < 0 + 1e-10;
}

static bool right_boundary(Mesh& Reader, int index_line) {
    // Returns TRUE if the LINE element at index index_line is on the right boundary
    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][0] > 1 - 1e-10\
        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][0] > 1 - 1e-10;
}

static bool bottom_boundary(Mesh& Reader, int index_line) {
    // Returns TRUE if the LINE element at index index_line is on the bottom boundary
    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][1] < 0 + 1e-10\
        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][1] < 0 + 1e-10;
}

static bool top_boundary(Mesh& Reader, int index_line) {
    // Returns TRUE if the LINE element at index index_line is on the top boundary
    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][1] > 1 - 1e-10\
        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][1] > 1 - 1e-10;
}

static bool all_boundary(Mesh& Reader, int index_node) {
    // Returns TRUE if the NODE at index index_node in on a boundary of the domain
    return (Reader.Nodes[index_node][0] < 0 + 1e-10 || Reader.Nodes[index_node][0] > 1 - 1e-10 || \
        Reader.Nodes[index_node][1] < 0 + 1e-10 || Reader.Nodes[index_node][1] > 1 - 1e-10);
}

static bool left_boundary_node(Mesh& Reader, int index_node) {
    // Returns TRUE if the NODE at index index_node in on a boundary of the domain
    return (Reader.Nodes[index_node][0] < 0 + 1e-10 || Reader.Nodes[index_node][1] < 0 + 1e-10);
}

//---------------------------------------------


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
    for (int i = 0; i < array_to_store.size(); i++) {
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

    Matrix_Builder MBuild;

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

    Volume_Inner_grad_Integral IntegrateK;
    MBuild.build_matrix(Reader, K, IntegrateK, K_function);

    cout << "Building force vector" << endl;
    Volume_Vector_Integral IntegrateF;
    MBuild.build_vector(Reader, F, IntegrateF, f_function);
    // Put F into previous F vector
    std::copy(F.begin(), F.end(), F_previous.begin());


    store_1d_vector_in_file("F", F);

    // Reset F vector before next time step
    std::fill(F.begin(), F.end(), 0.0);

    // For time dependent problems
    // Precomputing the constant matrices C and K
    cout << "Computing damping matrix" << endl;
    Volume_Matrix_Integral IntegrateC;
    MBuild.build_matrix(Reader, C, IntegrateC, C_function);

    store_2d_vector_in_file("C", C);
    store_2d_vector_in_file("K", K);

    // First, we calculate the initial condition at t=0
    MBuild.build_initial_T(Reader, dP_previous, T0);
    store_1d_vector_in_file("result_t_0", dP_previous);

    // // Then we loop on each timestep
    for (int i = 0; i < n_time; i++) {
        double t = (i + 1) * dt;
        cout << "Timestep t = " << t << endl;
        std::vector<std::vector<double>> A_matrix(n_dof, std::vector<double>(n_dof));
        std::vector<double> b_vector(n_dof);
        std::vector<double> params_for_fct = { t };

        MBuild.build_vector(Reader, F, IntegrateF, f_function);

        

        build_A_matrix(A_matrix);
        build_b_vector(b_vector);

        // Neumann boundary imposition
        cout << "Neumann BC ";
        Boundary_Vector_Integral Neumann_integrator;
        MBuild.neumann_BC(Reader, b_vector, Neumann_integrator, neumann_BC_fct, left_boundary, params_for_fct);
        store_1d_vector_in_file("F", F);

        // Robin boundary imposition
        cout << "Robin BC ";
        Boundary_Matrix_Integral Robin_integrator;
        MBuild.robin_BC(Reader, b_vector, A_matrix, Neumann_integrator, Robin_integrator, Robin_BC_vect_fct, Robin_BC_mat_fct, right_boundary, params_for_fct);

        // Dirichlet boundary imposition
        //MBuild.dirichlet_BC(Reader, A_matrix, b_vector, u0_fct, params_dirichlet, left_boundary_node);

        params_for_fct.clear();

        Solver.solve_system(A_matrix, b_vector, dP);
        store_1d_vector_in_file("b", b_vector);
        store_2d_vector_in_file("A", A_matrix);


        A_matrix.clear();
        b_vector.clear();

        store_1d_vector_in_file("result_t_" + to_string(t), dP);
        dP_previous.swap(dP);
        F_previous.swap(F);
        std::fill(F.begin(), F.end(), 0.0);
    }
    return 0;
};
