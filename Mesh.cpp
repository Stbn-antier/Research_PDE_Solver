#include "Mesh.h"

void Mesh::MeshReaderMSH(string filename)
{
    ifstream myfile;
    myfile.open(filename);

    string line;
    while (getline(myfile, line)) {
        if (line == "$MeshFormat") {
            getline(myfile, line);
            istringstream iss(line);
            iss >> file_version >> file_type >> data_size;
            getline(myfile, line);
        }
        if (line == "$Nodes") {
            getline(myfile, line);
            istringstream iss(line);
            int num_blocks;
            int min_tag;
            int max_tag;
            iss >> num_blocks >> num_nodes >> min_tag >> max_tag;
            for (int i = 0; i < num_blocks; i++) {
                getline(myfile, line);
                istringstream iss(line);
                int ent_dim;
                int ent_tag;
                bool param;
                int n_node_in_block;
                iss >> ent_dim >> ent_tag >> param >> n_node_in_block;

                for (int j = 0; j < n_node_in_block; j++) {
                    getline(myfile, line);
                }
                for (int j = 0; j < n_node_in_block; j++) {
                    getline(myfile, line);
                    istringstream iss(line);
                    double x; double y; double z;
                    iss >> x >> y >> z;
                    vector<double> node{ x, y, z };
                    Nodes.push_back(node);
                }
            }
        }
        if (line == "$Elements") {
            getline(myfile, line);
            istringstream iss(line);
            int num_blocks;
            int num_elems;
            int min_tag;
            int max_tag;
            iss >> num_blocks >> num_elems >> min_tag >> max_tag;

            for (int i = 0; i < num_blocks; i++) {
                getline(myfile, line);
                istringstream iss(line);
                int ent_dim;
                int ent_tag;
                int elem_type;
                int n_elem_in_block;
                iss >> ent_dim >> ent_tag >> elem_type >> n_elem_in_block;

                // Quad4 elements, 2D with 4 nodes
                if (elem_type == 3) {
                    num_Elems["quad"] += n_elem_in_block;
                    for (int j = 0; j < n_elem_in_block; j++) {
                        getline(myfile, line);
                        istringstream iss(line);
                        int id; int id_1; int id_2; int id_3; int id_4;
                        iss >> id >> id_1 >> id_2 >> id_3 >> id_4;
                        Element Elem(ent_tag, vector<int> {id_1 - 1, id_2 - 1, id_3 - 1, id_4 - 1});
                        Elems["quad"].push_back(Elem);
                    }
                }
                else if (elem_type == 1) {
                    num_Elems["line"] += n_elem_in_block;
                    for (int j = 0; j < n_elem_in_block; j++) {
                        getline(myfile, line);
                        istringstream iss(line);
                        int id; int id_1; int id_2;
                        iss >> id >> id_1 >> id_2;
                        Element Elem(ent_tag, vector<int> {id_1 - 1, id_2 - 1});
                        Elems["line"].push_back(Elem);
                    }
                }
                else {
                    for (int j = 0; j < n_elem_in_block; j++) {
                        getline(myfile, line);
                    }
                }
            }
        }
    }
}

int Mesh::FindBoundaryElementWithinMesh(Element elem, std::string searchtype) {
    // Size of the Higher order element, and lower order element
    int size_dplus_element = Elems[searchtype][0].Nodes.size();
    int size_search_element = elem.Nodes.size();
    for (int c = 0; c < num_Elems[searchtype]; c++) {
        int sim_counter = 0;
        for (int i = 0; i < size_dplus_element; i++) {
            for (int j = 0; j < size_search_element; j++) {
                if (Elems[searchtype][c].Nodes[i] == elem.Nodes[j]) {
                    sim_counter++;
                }
            }
        }
        if (sim_counter == size_search_element) {
            return c;
        }
    }
    throw(std::logic_error("The boundary isn't in the Elements"));
}

int Mesh::FindPositionOfBoundaryWithinElement(Element elem_boundary, Element elem_domain) {
    // For now, this function will only be able to find the boundary on a quad4 or tri3 element
    // where the boundary element is a line
    // Size of both elements
    int size_bound_element = elem_boundary.Nodes.size();
    int size_domain_element = elem_domain.Nodes.size();
    std::set<int> id_array;
    for (int i = 0; i < size_domain_element; i++) {
        for (int j = 0; j < size_bound_element; j++) {
            if (elem_domain.Nodes[i] == elem_boundary.Nodes[j]) {
                id_array.insert(i);
            }
        }
    }
    int bound_value;
    // For quad and tri element
    if (id_array == std::set<int> {0, 1}) {
        bound_value = 0;
    }
    else if (id_array == std::set<int> {1, 2}) {
        bound_value = 1;
    }
    else if (id_array == std::set<int> {2, 3}) {
        bound_value = 2;
    }
    // For quad element
    else if (id_array == std::set<int> {0, 3}) {
        bound_value = 3;
    }
    // For tri element
    else if (id_array == std::set<int> {0, 2}) {
        bound_value = 2;
    }
    id_array.clear();
    return bound_value;
}

Mesh::~Mesh()
{
}