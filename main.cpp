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

const double t_fin = 20; // Final time of simulation
const double n_time = 200; // Number of time steps not including t=0
const double dt = t_fin / n_time; // time step

std::vector<std::vector<double>> K; // Stiffness matrix
std::vector<double> F; // Force vector
std::vector<double> F_previous; // Force vector
std::vector<double> dP; // Solution vector
std::vector<double> dP_previous; // Solution vector at time t-1

std::vector<std::vector<double>> C; // Damping matrix for time stepping

//---------------------------------------------
// Constants for the different functions

const double k_iso = 44.5; // Heat transfer coefficient in W/(m.K)
const double rho = 7850; // Density of the material in kg/m^3
const double Cp = 475; // Heat capacity in J/(kg.K)

const double T_ext = 300; // Outside temperature in K
const double T_0 = 1100; // Initial temperature in K

const double h_conv = 730; // Convective heat coeff

//---------------------------------------------

const double Lxmin = -0.0125;
const double Lxmax = 0.0125;
const double Lymin = -0.0125;
const double Lymax = 0.0125;
const double Lzmin = 0;
const double Lzmax = 0.1;

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
    return 0;
}

//static double u0_fct(double x, double y, double z, std::vector<double>& params) {
//    return 0;
//}

// Boundary condition functions

static double Neumann_insulation(double x, double y, double z, std::vector<double>& params) {
    return 0;
}

static double Neumann_in_flux(double x, double y, double z, std::vector<double>& params) {
    //if ((-0.3 <= y && y <= 0.3) && (-0.3 <= z && z <= 0.3)) { // && (params[0] >= 7200)
    //    return q;
    //}
    //else {
    //    return 0;
    //}
    return q;
}

static double Robin_matrix(double x, double y, double z, std::vector<double>& params) {
    return h;
}

static double Robin_vector(double x, double y, double z, std::vector<double>& params) {
    return h * T_ext;
}

// Initial time function

static double T0(std::vector<double>& coordinates) {
    return T_0;
}


//---------------------------------------------
// Functions defining where the boundaries are

// Boundary for 2D elements
static bool left_surf(Mesh& Reader, int index_quad) {
    // Returns TRUE if the quad is on the left bound
    
    bool res = 1;
    for (int i = 0; i < 4; i++) { //Iterate on nodes of quad element
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][0] < Lxmin + 1e-3);
    }
    return res;
}

static bool right_surf(Mesh& Reader, int index_quad) {
    bool res = 1;
    for (int i = 0; i < 4; i++) { //Iterate on nodes of quad element
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][0] > Lxmax - 1e-3);
    }
    return res;
}

static bool top_surf(Mesh& Reader, int index_quad) {
    bool res = 1;
    for (int i = 0; i < 4; i++) { //Iterate on nodes of quad element
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][2] > Lzmax - 1e-3);
    }
    return res;
}

static bool bottom_surf(Mesh& Reader, int index_quad) {
    bool res = 1;
    for (int i = 0; i < 4; i++) { //Iterate on nodes of quad element
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][2] < Lzmin + 1e-3);
    }
    return res;
}

static bool front_surf(Mesh& Reader, int index_quad) {
    bool res = 1;
    for (int i = 0; i < 4; i++) { //Iterate on nodes of quad element
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][1] > Lymax - 1e-3);
    }
    return res;
}

static bool back_surf(Mesh& Reader, int index_quad) {
    bool res = 1;
    for (int i = 0; i < 4; i++) { //Iterate on nodes of quad element
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][1] < Lymin + 1e-3);
    }
    return res;
}

static bool insulation_bound(Mesh& Reader, int index_quad) {
    return (top_surf(Reader, index_quad) || bottom_surf(Reader, index_quad) || front_surf(Reader, index_quad) || back_surf(Reader, index_quad));
}

// Boundary for Nodes
static bool all_boundary(Mesh& Reader, int index_node) {
    // Returns TRUE if the NODE at index index_node in on a boundary of the domain
    return (\
        Reader.Nodes[index_node][0] < Lxmin + 1e-3 || Reader.Nodes[index_node][0] > Lxmax - 1e-3 || \
        Reader.Nodes[index_node][1] < Lymin + 1e-3 || Reader.Nodes[index_node][1] > Lymax - 1e-3 || \
        Reader.Nodes[index_node][2] < Lzmin + 1e-3 || Reader.Nodes[index_node][2] > Lzmax - 1e-3);
}

static bool bottom(Mesh& Reader, int index_node) {
    return (Reader.Nodes[index_node][2] < Lzmin + 1e-3);
}

static bool bottom_top(Mesh& Reader, int index_node) {
    return (Reader.Nodes[index_node][2] < Lzmin + 1e-3 || Reader.Nodes[index_node][2] > Lzmax - 1e-3);
}

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

void store_2d_vector_in_file(std::string filename, vector<vector<double>>& array_to_store, int lines_to_take) {
    ofstream myfile;
    myfile.open(filename + ".txt");
    for (int i = 0; i < lines_to_take; i++) {
        for (int j = 0; j < lines_to_take; j++) {
            myfile << array_to_store[i][j] << ' ';
        }
        myfile << endl;
    }
    myfile.close();
}

int main() {

    // Read mesh in file, prepare result writer
    Mesh Reader;
    Reader.MeshReaderMSH("Mesh/mesh_10x10x10.msh");

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



    // Building constant matrixes and vectors
    cout << "Building stiffness matrix" << endl;
    Volume_Inner_grad_Integral3D IntegrateK;
    MBuild.build_matrix(Reader, K, IntegrateK, K_function);

    cout << "Building force vector" << endl;
    Volume_Vector_Integral3D IntegrateF;
    MBuild.build_vector(Reader, F, IntegrateF, f_function);

    // Neumann BC
    //std::cout << "Neumann BC " << std::endl;
    //std::vector<double> params_neumann = { 0.0 };
    //Boundary_Vector_Integral3D Neumann_integrator;
    //MBuild.neumann_BC(Reader, F, Neumann_integrator, Neumann_insulation, insulation_bound, params_neumann);
    //MBuild.neumann_BC(Reader, F, Neumann_integrator, Neumann_in_flux, left_surf, params_neumann);

    //// Robin BC
    //std::cout << "Robin BC " << std::endl;
    //std::vector<double> params_robin = { 0.0 };
    //Boundary_Matrix_Integral3D Robin_integrator;
    //MBuild.robin_BC(Reader, F, K, Neumann_integrator, Robin_integrator, \
    //    Robin_vector, Robin_matrix, right_surf, params_robin);

    
    // Put F into previous F vector
    std::copy(F.begin(), F.end(), F_previous.begin());
    //store_1d_vector_in_file("results/F", F);

    // Reset F vector before next time step
    std::fill(F.begin(), F.end(), 0.0);

    // For time dependent problems
    // Precomputing the constant matrices C and K
    cout << "Computing damping matrix" << endl;
    Volume_Matrix_Integral3D IntegrateC;
    MBuild.build_matrix(Reader, C, IntegrateC, C_function);


    store_2d_vector_in_file("results/C", C, 10);
    store_2d_vector_in_file("results/K", K, 10);


    //// Dirichlet boundary imposition
    //std::vector<double> params_dirichlet = { 0.0 };
    ////MBuild.dirichlet_BC(Reader, K, F_previous, u0_fct, params_dirichlet, all_boundary);
    //MBuild.dirichlet_BC(Reader, K, F_previous, u0_fct, params_dirichlet, bottom_top);


    //Solver.solve_system(K, F_previous, dP);

    //store_1d_vector_in_file("results/T", dP);
    //store_1d_vector_in_binary_file("results/T", dP);
    //Writer.Write_time_step(Reader, 1.0, "T", dP);

    // First, we calculate the initial condition at t=0
    MBuild.build_initial_T(Reader, dP_previous, T0);
    Writer.Write_time_step(Reader, 0, "T", dP_previous);
    store_1d_vector_in_binary_file("results/T0", dP_previous);


    for (int i = 0; i < n_time; i++) {
        double t = (i + 1) * dt; // Time in seconds
        std::cout << "Timestep t = " << t << " s" << endl;

        std::vector<std::vector<double>> A_matrix(n_dof, std::vector<double>(n_dof));
        std::vector<double> b_vector(n_dof);
        
        std::cout << "Build force vector" << std::endl;
        MBuild.build_vector(Reader, F, IntegrateF, f_function);
        
        // Because K and C are constant
        build_A_matrix(A_matrix);
        build_b_vector(b_vector);

        // Neumann BC
        std::cout << "Neumann BC " << std::endl;
        std::vector<double> params_neumann = { t };
        Boundary_Vector_Integral3D Neumann_integrator;
        MBuild.neumann_BC(Reader, b_vector, Neumann_integrator, Neumann_insulation, insulation_bound, params_neumann);
        MBuild.neumann_BC(Reader, b_vector, Neumann_integrator, Neumann_in_flux, left_surf, params_neumann);

        // Robin BC
        std::cout << "Robin BC " << std::endl;
        std::vector<double> params_robin = { 0.0 };
        Boundary_Matrix_Integral3D Robin_integrator;
        MBuild.robin_BC(Reader, b_vector, A_matrix, Neumann_integrator, Robin_integrator, \
            Robin_vector, Robin_matrix, right_surf, params_robin);

        std::cout << "Solving system" << std::endl;
        Solver.solve_system_Cholesky(A_matrix, b_vector, dP);


        A_matrix.clear();
        b_vector.clear();

        Writer.Write_time_step(Reader, t, "T", dP);
        store_1d_vector_in_binary_file("results/T"+std::to_string(i+1), dP_previous);


        dP_previous.swap(dP);
        F_previous.swap(F);
        std::fill(F.begin(), F.end(), 0.0);
    }

    return 0;
};
