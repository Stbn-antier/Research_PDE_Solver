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
#include "Matrix.h"

using namespace std;

const double t_fin = 10; // Final time of simulation
const int n_time = 50; // Number of time steps not including t=0
const double dt = t_fin / n_time; // time step

//std::vector<std::vector<double>> K; // Stiffness matrix
//std::vector<double> F; // Force vector
//std::vector<double> F_previous; // Force vector
//std::vector<double> dP; // Solution vector
//std::vector<double> dP_previous; // Solution vector at time t-1
//
//std::vector<std::vector<double>> C; // Damping matrix for time stepping

//---------------------------------------------
// Constants for the different functions

const double k_iso = 44.5; // Heat transfer coefficient in W/(m.K)
const double rho = 7850; // Density of the material in kg/m^3
const double Cp = 475; // Heat capacity in J/(kg.K)

const double T_ext = 300; // Outside temperature in K
const double T_0 = 1100; // Initial temperature in K

const double hconv = 730; // Convective heat coeff
//const double hnb = 15000; // Nucleation boiling heat transfer
const double hfb = 100; // Heat transfer with oil vapor

double hnb_list[3] = { 10000 , 15000 , 20000 };
double Lnb_list[3] = { 0.005 , 0.01 , 0.015 };
double Vconv_list[3] = { 0.02 , 0.03 , 0.04 }; // Speed of the nucleation bubbling layer displacement

const double stand_dev = 0.005; // Standard deviation of the gaussian distribution used
const int num_loops = 100;

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

static double h(double x, double y, double z, std::vector<double>& params) {
    // Params vector stores information for planes at the boundary
    // params[0] -> n1x, params[1] -> n1y, params[2] -> d1, params[3] -> n2x, params[4] -> n2y, params[5] -> d2
    // params[6] -> hnb
    double above1 = params[0] * x + params[1] * y + z - params[2];
    double above2 = params[3] * x + params[4] * y + z - params[5];
    // If above1>0, we are "above" plane 1, same for above2>0 and plane 2

    if (above1 <= 0 && above2 <= 0) {
        return hconv;
    }
    else if (above1 > 0 && above2 > 0) {
        return hfb;
    }
    else {
        return params[6];
    }
}

// Boundary condition functions

static double Neumann_insulation(double x, double y, double z, std::vector<double>& params) {
    return 0;
}

static double Robin_matrix(double x, double y, double z, std::vector<double>& params) {
    return h(x,y,z,params);
}

static double Robin_vector(double x, double y, double z, std::vector<double>& params) {
    return h(x, y, z, params) * T_ext;
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
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][0] < Lxmin + 1e-10);
    }
    return res;
}

static bool right_surf(Mesh& Reader, int index_quad) {
    bool res = 1;
    for (int i = 0; i < 4; i++) { //Iterate on nodes of quad element
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][0] > Lxmax - 1e-10);
    }
    return res;
}

static bool top_surf(Mesh& Reader, int index_quad) {
    bool res = 1;
    for (int i = 0; i < 4; i++) { //Iterate on nodes of quad element
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][2] > Lzmax - 1e-10);
    }
    return res;
}

static bool bottom_surf(Mesh& Reader, int index_quad) {
    bool res = 1;
    for (int i = 0; i < 4; i++) { //Iterate on nodes of quad element
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][2] < Lzmin + 1e-10);
    }
    return res;
}

static bool front_surf(Mesh& Reader, int index_quad) {
    bool res = 1;
    for (int i = 0; i < 4; i++) { //Iterate on nodes of quad element
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][1] > Lymax - 1e-10);
    }
    return res;
}

static bool back_surf(Mesh& Reader, int index_quad) {
    bool res = 1;
    for (int i = 0; i < 4; i++) { //Iterate on nodes of quad element
        res &= (Reader.Nodes[Reader.Elems["quad"][index_quad].Nodes[i]][1] < Lymin + 1e-10);
    }
    return res;
}

static bool insulation_bound(Mesh& Reader, int index_quad) {
    return (top_surf(Reader, index_quad) || bottom_surf(Reader, index_quad) || front_surf(Reader, index_quad) || back_surf(Reader, index_quad));
}

static bool all_surf(Mesh& Reader, int index_quad) {
    return (top_surf(Reader, index_quad) || bottom_surf(Reader, index_quad)\
        || front_surf(Reader, index_quad) || back_surf(Reader, index_quad)\
        || left_surf(Reader, index_quad) || right_surf(Reader, index_quad));
}

// Boundary for Nodes
static bool all_boundary(Mesh& Reader, int index_node) {
    // Returns TRUE if the NODE at index index_node in on a boundary of the domain
    return (\
        Reader.Nodes[index_node][0] < Lxmin + 1e-10 || Reader.Nodes[index_node][0] > Lxmax - 1e-10 || \
        Reader.Nodes[index_node][1] < Lymin + 1e-10 || Reader.Nodes[index_node][1] > Lymax - 1e-10 || \
        Reader.Nodes[index_node][2] < Lzmin + 1e-10 || Reader.Nodes[index_node][2] > Lzmax - 1e-10);
}

static bool bottom(Mesh& Reader, int index_node) {
    return (Reader.Nodes[index_node][2] < Lzmin + 1e-10);
}

static bool bottom_top(Mesh& Reader, int index_node) {
    return (Reader.Nodes[index_node][2] < Lzmin + 1e-10 || Reader.Nodes[index_node][2] > Lzmax - 1e-10);
}

//---------------------------------------------

// Functions for storing data in files
void store_1d_vector_in_file(std::string filename, vector<double>& array_to_store) {
    ofstream myfile;
    myfile.open(filename + ".txt");
    for (int i = 0; i < array_to_store.size(); i++) {
        myfile << array_to_store[i] << endl;
    }
    myfile.close();
}

void store_1d_vector_in_file(std::string filename, DataVector& array_to_store) {
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

void store_1d_vector_in_binary_file(std::string filename, DataVector& array_to_store) {
    std::ofstream myfile(filename + ".bin", std::ios::binary);
    myfile.write(reinterpret_cast<const char*>(array_to_store.data_vec()), array_to_store.size() * sizeof(double));

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

//---------------------------------------------------------------

int main() {
    // RNG setup
    std::random_device rd{};
    std::mt19937 gen{ rd() };

    std::normal_distribution<double> gauss{ 0, stand_dev };
    auto random_float = [&gen, &gauss] {return gauss(gen); };
    std::vector<std::vector<double>> random_process(num_loops, std::vector<double>(n_time));
    for (int i = 0; i < n_time; i++) {
        for (int j = 0; j < num_loops; j++) {
            random_process[j][i] = random_float();
        }
    }
    store_2d_vector_in_file("results/random_process", random_process);

    // Read mesh in file, prepare result writer
    Mesh Reader;
    Reader.MeshReaderMSH("Mesh/mesh_3d_quench.msh");

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
    
    // Create vectors and matrices to hold data
    int n_dof = Reader.num_nodes;
    COOMatrix K(true);
    COOMatrix C(true);

    // Building stiffness and damping matrices
    cout << "Building stiffness matrix" << endl;
    Volume_Inner_grad_Integral3D IntegrateK;
    MBuild.build_matrix(Reader, K, IntegrateK, K_function);

    cout << "Computing damping matrix" << endl;
    Volume_Matrix_Integral3D IntegrateC;
    MBuild.build_matrix(Reader, C, IntegrateC, C_function);

    Volume_Vector_Integral3D IntegrateF;

    // Make C and K into CSR Matrices
    CSRMatrix K_CSR(K);
    CSRMatrix C_CSR(C);


    for (int i_h = 0; i_h < 3; i_h++) {
        for (int i_L = 0; i_L < 3; i_L++) {
            for (int i_Ldot = 0; i_Ldot < 3; i_Ldot++) {
                for (int i_loop = 0; i_loop < num_loops; i_loop++) {
                    std::string path_storage = "results/loop_" + to_string(i_h * 9 * num_loops + i_L * 3 * num_loops + i_Ldot * num_loops + i_loop + 1) + "/";
                    std::filesystem::create_directory(path_storage);

                    // We set the parameters for now:
                    double hnb = hnb_list[i_h];
                    double Lnb = Lnb_list[i_L];
                    double Vconv = Vconv_list[i_Ldot];

                    // Vectors for solving
                    DataVector F(n_dof);
                    DataVector F_previous(n_dof);
                    DataVector dP(n_dof);
                    DataVector dP_previous(n_dof);

                    std::vector<double> T_C(n_time + 1); // Storing the values of temperature at center
                    //std::cout << Reader.Nodes[3381][0] << " " << Reader.Nodes[3381][1] << " " << Reader.Nodes[3381][2] << std::endl;
                    // for 10x10x40 mesh Center of sample is node 3381

                    // Build force vector for timestep 0
                    MBuild.build_vector(Reader, F, IntegrateF, f_function);
                    // Put F into previous F vector
                    std::copy(F.begin(), F.end(), F_previous.begin());
                    // Reset F vector before next time step
                    std::fill(F.begin(), F.end(), 0.0);


                    // First, we calculate the initial condition at t=0
                    MBuild.build_initial_T(Reader, dP_previous, T0);
                    T_C[0] = dP_previous[3381];
                    store_1d_vector_in_binary_file(path_storage + "T0", dP_previous);

                    std::vector<double> Lconv_vector(n_time + 1); // Storage for the realizations of random

                    for (int i = 0; i < n_time; i++) {
                        double t = (i + 1) * dt;
                        double gaussian_number = random_process[i_loop][i];
                        Lconv_vector[i + 1] = Lconv_vector[i] + dt * Vconv + gaussian_number;

                        cout << "Loop number n=" << i_loop + 1 << "/" << num_loops << ", Timestep t = " << t << endl;
                        cout << "Params : hnb=" << hnb << " Lnb=" << Lnb << " Vconv=" << Vconv << endl;
                        cout << "Random number : " << gaussian_number << endl;

                        std::cout << "Build force vector" << std::endl;
                        MBuild.build_vector(Reader, F, IntegrateF, f_function);

                        // Build A_matrix and b_vector
                        COOMatrix A_matrix(C_CSR + dt / 2 * K_CSR);
                        DataVector b_vector = (C_CSR - dt / 2 * K_CSR) * dP_previous + dt / 2 * (F + F_previous);


                        // Robin BC
                        std::cout << "Robin BC " << std::endl;
                        // Params vector stores information for planes at the boundary
                        // params[0] -> n1x, params[1] -> n1y, params[2] -> d1, params[3] -> n2x, params[4] -> n2y, params[5] -> d2, params[6] -> hnb
                        std::vector<double> params_robin = { 0, 0, Lconv_vector[i + 1], 0, 0 , Lconv_vector[i + 1] + Lnb, hnb };
                        Boundary_Vector_Integral3D Neumann_integrator;
                        Boundary_Matrix_Integral3D Robin_integrator;
                        MBuild.robin_BC(Reader, b_vector, A_matrix, Neumann_integrator, Robin_integrator, \
                            Robin_vector, Robin_matrix, all_surf, params_robin);

                        std::cout << "Solving system" << std::endl;
                        // First convert from COO to CSR for solver
                        CSRMatrix A_CSR(A_matrix);
                        // Solve system with iterative PCCG solver
                        dP = Solver.PCCG(A_CSR, b_vector, dP_previous);


                        // Store solution
                        store_1d_vector_in_binary_file(path_storage + "T" + std::to_string(i + 1), dP);
                        T_C[i + 1] = dP[3381];


                        dP_previous.swap(dP);
                        F_previous.swap(F);
                        std::fill(F.begin(), F.end(), 0.0);
                    }
                    store_1d_vector_in_file(path_storage + "Lconv", Lconv_vector);
                    store_1d_vector_in_file(path_storage + "T_C", T_C);
                }
            }
        }
    }

    return 0;
};
