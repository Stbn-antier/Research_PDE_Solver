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
#include "Shape_fct_1D.h"
#include "Solve_matrix_system.h"
#include "Integrals.h"
#include "Matrix_Assembly.h"

using namespace std;

int n_elem = 464;
const int n_sf = 4; // Number of shape function per element
const int n_dof = 1411; // Number of nodes | à modifier à chaque fois, 505 pour t2, 143 pour t1

const int report_step = 4; // Report progress at every 1/report_step % of progress

const double t_fin = 10; // Final time of simulation
const double n_time = 50; // Number of time steps not including t=0
const double dt = t_fin / n_time; // time step

std::vector<std::vector<double>> K(n_dof, std::vector<double>(n_dof)); // Stiffness matrix
std::vector<double> F(n_dof); // Force vector
std::vector<double> F_previous(n_dof); // Force vector
std::vector<double> dP(n_dof); // Solution vector
std::vector<double> dP_previous(n_dof); // Solution vector at time t-1

std::vector<std::vector<double>> C(n_dof, std::vector<double>(n_dof)); // Damping matrix for time stepping

//---------------------------------------------
// Constants for the different functions
const double lx2 = 0.0125; // Width of the sample (half of total width)
const double Ly = 0.1; // Height of the sample

const double k_iso = 44.5; // Heat transfer coefficient in W/m²
const double rho = 7850; // Density of the material in kg/m^3
const double Cp = 475; // Heat capacity in J/(kg.K)

const double T_ext = 300; // Outside temperature in K
const double T_0 = 1100; // Initial temperature in K

const double hconv = 730;
//const double hnb = 15000;
const double hfb = 140;

const double Lnb = 0.03;
//const double ynb = 0.02;
const double xfb = 0.005;

//---------------------------------------------
// Parameters to vary
//const double tau = 8; // Time for the heat transfer to go throughout the sample
std::vector<double> tau_list = { 2 , 4 , 6 , 8 };
std::vector<double> hnb_list = { 10000 , 12500 , 15000 , 17500 , 20000 };
std::vector<double> Lnb_list = { 0.005 , 0.01, 0.02 , 0.03 , 0.04 };
const double stand_dev = 5; // Standard deviation of the gaussian distribution used


const double h_coef_right(double y, std::vector<double>& params) {
    // Heat transfer coefficient for the right boundary in W/(m².K)
    // Additional parameters are :
    // params[0] -> t, params[1] -> random_float, params[2] -> tau
    //

    double ynb = Ly * (params[0]) / (params[2]);

    if (y < ynb) {
        return hconv;
    }
    else if (ynb <= y && y <= ynb + params[4]) {
        return params[3];
    }
    else if (ynb + params[4] < y) {
        return hfb;
    }
    return NAN;
}

const double h_coef_top(double x, std::vector<double>& params) {
    // Heat transfer coefficient for the right boundary in W/(m².K)
    if (x <= xfb) {
        return hconv;
    }
    else if (xfb < x) {
        return hfb;
    }
    return NAN;
}

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

// Boundary condition functions

static double Robin_BC_vect_fct_right(double x, double y, std::vector<double>& params) {
    return h_coef_right(y, params) * T_ext;
}

static double Robin_BC_mat_fct_right(double x, double y, std::vector<double>& params) {
    return h_coef_right(y, params);
}

static double Robin_BC_vect_fct_top(double x, double y, std::vector<double>& params) {
    return h_coef_top(x, params) * T_ext;
}

static double Robin_BC_mat_fct_top(double x, double y, std::vector<double>& params) {
    return h_coef_top(x, params);
}

static double Robin_BC_vect_fct_bottom(double x, double y, std::vector<double>& params) {
    return hconv * T_ext;
}

static double Robin_BC_mat_fct_bottom(double x, double y, std::vector<double>& params) {
    return hconv;
}

// Initial time function

static double T0(std::vector<double>& coordinates) {
    return T_0;
}


//---------------------------------------------
// Functions defining where the boundaries are

static bool left_boundary(Mesh& Reader, int index_line) {
    // Returns TRUE if the LINE element at index index_line is on the left boundary
    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][0] < 0 + 1e-10\
        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][0] < 0 + 1e-10;
}

static bool right_boundary(Mesh& Reader, int index_line) {
    // Returns TRUE if the LINE element at index index_line is on the right boundary
    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][0] > lx2 - 1e-10\
        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][0] > lx2 - 1e-10;
}

static bool bottom_boundary(Mesh& Reader, int index_line) {
    // Returns TRUE if the LINE element at index index_line is on the bottom boundary
    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][1] < 0 + 1e-10\
        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][1] < 0 + 1e-10;
}

static bool top_boundary(Mesh& Reader, int index_line) {
    // Returns TRUE if the LINE element at index index_line is on the top boundary
    return Reader.Nodes[Reader.Elems["line"][index_line].Nodes[0]][1] > Ly - 1e-10\
        && Reader.Nodes[Reader.Elems["line"][index_line].Nodes[1]][1] > Ly - 1e-10;
}

static bool all_boundary(Mesh& Reader, int index_node) {
    // Returns TRUE if the NODE at index index_node in on a boundary of the domain
    return (Reader.Nodes[index_node][0] < 0 + 1e-10 || Reader.Nodes[index_node][0] > lx2 - 1e-10 || \
        Reader.Nodes[index_node][1] < 0 + 1e-10 || Reader.Nodes[index_node][1] > Ly - 1e-10);
}

//static bool left_boundary_node(Mesh& Reader, int index_node) {
//    // Returns TRUE if the NODE at index index_node in on a boundary of the domain
//    return (Reader.Nodes[index_node][0] < 0 + 1e-10 || Reader.Nodes[index_node][1] < 0 + 1e-10);
//}

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
    myfile.open(filename + ".txt");
    for (int i = 0; i < array_to_store.size(); i++) {
        myfile << array_to_store[i] << endl;
    }
    myfile.close();
}

void store_2d_vector_in_file(std::string filename, vector<vector<double>>& array_to_store) {
    ofstream myfile;
    myfile.open(filename + ".txt");
    for (int i = 0; i < n_dof; i++) {
        for (int j = 0; j < n_dof; j++) {
            myfile << array_to_store[i][j] << ' ';
        }
        myfile << endl;
    }
    myfile.close();
}

int main() {
    // RNG setup
    std::random_device rd{};
    std::mt19937 gen{ rd() };

    std::normal_distribution<double> gauss{0, stand_dev};
    auto random_float = [&gen, &gauss] {return gauss(gen);};

    // Mesh Setup, matrix Builder Setup
    Mesh Reader;
    Reader.MeshReaderMSH("mesh_rectangle.msh");

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


    // For plotting values at specific points, first setup vectors and find the DOFs of specific points
    int dof_a = Reader.FindDofFromCoords(0, 0.75 * Ly, 1e-6);
    cout << "Point A is DOF number " << dof_a << std::endl;
    int dof_b = Reader.FindDofFromCoords(lx2 / 2, 0.75 * Ly, 1e-6);
    cout << "Point B is DOF number " << dof_b << std::endl;
    int dof_c = Reader.FindDofFromCoords(0, Ly/2, 1e-6);
    cout << "Center of domain (Point C) is DOF number " << dof_c << std::endl;
    int dof_d = Reader.FindDofFromCoords(lx2 / 2, Ly/2, 1e-6);
    cout << "Point D is DOF number " << dof_d << std::endl;
    int dof_e = Reader.FindDofFromCoords(0, 0.25 * Ly, 1e-6);
    cout << "Point E is DOF number " << dof_e << std::endl;
    int dof_f = Reader.FindDofFromCoords(lx2 / 2, 0.25 * Ly, 1e-6);
    cout << "Point F is DOF number " << dof_f << std::endl;
    
    std::vector<double> probe_Ta(n_time + 1);
    std::vector<double> probe_Tb(n_time + 1);
    std::vector<double> probe_Tc(n_time + 1);
    std::vector<double> probe_Td(n_time + 1);
    std::vector<double> probe_Te(n_time + 1);
    std::vector<double> probe_Tf(n_time + 1);

    std::vector<double> rnd_vector(n_time + 1); // Storage for the realizations of random


    // Building constant matrixes and vectors
    cout << "Building stiffness matrix" << endl;
    Volume_Inner_grad_Integral IntegrateK;
    MBuild.build_matrix(Reader, K, IntegrateK, K_function);

    cout << "Building force vector" << endl;
    Volume_Vector_Integral IntegrateF;
    MBuild.build_vector(Reader, F, IntegrateF, f_function);
    // Put F into previous F vector
    std::copy(F.begin(), F.end(), F_previous.begin());


    //store_1d_vector_in_file("F", F);

    // Reset F vector before next time step
    std::fill(F.begin(), F.end(), 0.0);

    // For time dependent problems
    // Precomputing the constant matrices C and K
    cout << "Computing damping matrix" << endl;
    Volume_Matrix_Integral IntegrateC;
    MBuild.build_matrix(Reader, C, IntegrateC, C_function);

    //store_2d_vector_in_file("C", C);
    //store_2d_vector_in_file("K", K);

    // Looping on the different parameters used
    int total_length = tau_list.size() * hnb_list.size() * Lnb_list.size();
    for (int i_tau = 0; i_tau < tau_list.size(); i_tau++) {
        for (int i_hnb = 0; i_hnb < hnb_list.size(); i_hnb++) {
            for (int i_Lnb = 0; i_Lnb < Lnb_list.size(); i_Lnb++) {
                double tau = tau_list[i_tau];
                double hnb = hnb_list[i_hnb];
                double Lnb = Lnb_list[i_Lnb];

                std::string path_storage = "results/tau_" + std::format("{:.0f}", tau) + "_hnb_" + std::format("{:.0f}", hnb) + "_Lnb_" + std::format("{:.3f}", Lnb) + "/";
                std::filesystem::create_directory(path_storage);

                // First, we calculate the initial condition at t=0
                MBuild.build_initial_T(Reader, dP_previous, T0);
                store_1d_vector_in_file(path_storage + "result_step_0", dP_previous);
                probe_Ta[0] = dP_previous[dof_a];
                probe_Tb[0] = dP_previous[dof_b];
                probe_Tc[0] = dP_previous[dof_c];
                probe_Td[0] = dP_previous[dof_d];
                probe_Te[0] = dP_previous[dof_e];
                probe_Tf[0] = dP_previous[dof_f];
                rnd_vector[0] = 0;

                // // Then we loop on each timestep
                for (int i = 0; i < n_time; i++) {
                    double t = (i + 1) * dt;
                    double additional_param = random_float();
                    rnd_vector[i + 1] = additional_param;
                    cout << "Loop number " << i_Lnb+i_hnb+i_tau << "/" << total_length << " Timestep t = " << t << endl;
                    cout << "Random number : " << additional_param << endl;
                    std::vector<std::vector<double>> A_matrix(n_dof, std::vector<double>(n_dof));
                    std::vector<double> b_vector(n_dof);
                    std::vector<double> params_for_fct = { t , additional_param , tau , hnb , Lnb };

                    MBuild.build_vector(Reader, F, IntegrateF, f_function);



                    build_A_matrix(A_matrix);
                    build_b_vector(b_vector);

                    // Neumann boundary imposition
                    cout << "Neumann BC ";
                    Boundary_Vector_Integral Neumann_integrator;
                    //MBuild.neumann_BC(Reader, b_vector, Neumann_integrator, neumann_BC_fct, left_boundary, params_for_fct);
                    //store_1d_vector_in_file("F", F);

                    // Robin boundary imposition
                    cout << "Robin BC ";
                    Boundary_Matrix_Integral Robin_integrator;
                    // Top boundary
                    MBuild.robin_BC(Reader, b_vector, A_matrix, Neumann_integrator, Robin_integrator, \
                        Robin_BC_vect_fct_top, Robin_BC_mat_fct_top, top_boundary, params_for_fct);
                    // Right boundary
                    MBuild.robin_BC(Reader, b_vector, A_matrix, Neumann_integrator, Robin_integrator, \
                        Robin_BC_vect_fct_right, Robin_BC_mat_fct_right, right_boundary, params_for_fct);
                    // Bottom boundary
                    MBuild.robin_BC(Reader, b_vector, A_matrix, Neumann_integrator, Robin_integrator, \
                        Robin_BC_vect_fct_bottom, Robin_BC_mat_fct_bottom, bottom_boundary, params_for_fct);


                    // Dirichlet boundary imposition
                    //MBuild.dirichlet_BC(Reader, A_matrix, b_vector, u0_fct, params_dirichlet, left_boundary_node);

                    params_for_fct.clear();

                    Solver.solve_system(A_matrix, b_vector, dP);
                    //store_1d_vector_in_file("b", b_vector);
                    //store_2d_vector_in_file("A", A_matrix);

                    probe_Ta[i + 1] = dP[dof_a];
                    probe_Tb[i + 1] = dP[dof_b];
                    probe_Tc[i + 1] = dP[dof_c];
                    probe_Td[i + 1] = dP[dof_d];
                    probe_Te[i + 1] = dP[dof_e];
                    probe_Tf[i + 1] = dP[dof_f];


                    A_matrix.clear();
                    b_vector.clear();

                    store_1d_vector_in_file(path_storage + "/result_step_" + to_string(i + 1), dP);
                    dP_previous.swap(dP);
                    F_previous.swap(F);
                    std::fill(F.begin(), F.end(), 0.0);
                }

                store_1d_vector_in_file(path_storage + "T_A", probe_Ta);
                store_1d_vector_in_file(path_storage + "T_B", probe_Tb);
                store_1d_vector_in_file(path_storage + "T_C", probe_Tc);
                store_1d_vector_in_file(path_storage + "T_D", probe_Td);
                store_1d_vector_in_file(path_storage + "T_E", probe_Te);
                store_1d_vector_in_file(path_storage + "T_F", probe_Tf);

                //store_1d_vector_in_file(path_storage + "rnd_vector", rnd_vector);
            }
        }
    }
    return 0;
};
