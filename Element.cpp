#include "Element.h"

Element::Element(int tag, std::vector<int> node_id)
{
    elem_tag = tag;
    Nodes = node_id;
}

Element::~Element()
{
}

int Element::get_node(int i) {
    return Nodes[i];
}
