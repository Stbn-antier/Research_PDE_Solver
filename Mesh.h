#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <set>
#include "Element.h"

using namespace std;

class Mesh
{
private:
public:
    map<std::string, vector<Element>> Elems;
    map<std::string, int> num_Elems;
    vector<vector<double>> Nodes;
    double file_version = 0.;
    bool file_type = 0;
    int data_size = 0;
    int num_nodes = 0;
    int num_quad = 0;
    int shapefct_per_node = 4;
    void MeshReaderMSH(string filename);
    int FindBoundaryElementWithinMesh(Element elem, std::string search_type);
    int FindPositionOfBoundaryWithinElement(Element elem_boundary, Element elem_domain);
    int FindDofFromCoords(double x, double y, double epsilon);
    ~Mesh();
};

