#pragma once
#include <vector>

class Element
{
private:
public:
    int elem_tag;
    std::vector<int> Nodes;
    Element(int tag, std::vector<int> node_id);
    ~Element();
    int get_node(int i);
};

