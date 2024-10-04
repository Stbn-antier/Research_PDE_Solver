#include "Writer.h"

XML_Writer::XML_Writer(std::string path_storage, std::string filename)
{
	// Open file and write Header
	xmlfile.open(path_storage + filename + ".xdmf");
	xmlfile << "<?xml version=\"1.0\"?>" << std::endl;
	xmlfile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
	xmlfile << "<Xdmf Version=\"3.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">" << std::endl;
	xmlfile << "<Domain>" << std::endl;
	xmlfile << "<Grid Name=\"" + filename + "\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
}

XML_Writer::~XML_Writer()
{
	// Write footer and close file
	xmlfile << "</Grid>" << std::endl;
	xmlfile << "</Domain>" << std::endl;
	xmlfile << "</Xdmf>" << std::endl;
	xmlfile.close();
}

void XML_Writer::Write_topology_data(Mesh& Reader)
{
	// Write topology dataset (connectivity between nodes, definition of elements)
	xmlfile << "<Topology NumberOfElements=\"" + std::to_string(Reader.num_Elems["hexa"]) + "\" TopologyType=\"Hexahedron\" NodesPerElement=\"8\">" << std::endl;
	xmlfile << "<DataItem Dimensions=\"" + std::to_string(Reader.num_Elems["hexa"]) + " 8\" NumberType=\"UInt\" Format=\"XML\">" << std::endl;
	// Write data
	for (int i = 0; i < Reader.num_Elems["hexa"]; i++) {
		for (int j = 0; j < 8; j++) {
			xmlfile << std::to_string(Reader.Elems["hexa"][i].Nodes[j]) << " ";
		}
		xmlfile << std::endl;
	}

	xmlfile << "</DataItem>" << std::endl;
	xmlfile << "</Topology>" << std::endl;
}

void XML_Writer::Write_geometry_data(Mesh& Reader)
{
	xmlfile << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	xmlfile << "<DataItem Dimensions=\""+ std::to_string(Reader.num_nodes) + " 3\" Format=\"XML\">" << std::endl;
	// Write data
	for (int i = 0; i < Reader.num_nodes; i++) {
		for (int j = 0; j < 3; j++) {
			xmlfile << std::to_string(Reader.Nodes[i][j]) << " ";
		}
		xmlfile << std::endl;
	}
		
	xmlfile << "</DataItem>" << std::endl;
	xmlfile << "</Geometry>" << std::endl;
}

void XML_Writer::Write_time_step(Mesh& Reader, double time, std::string data_name, std::vector<double>& data)
{
	// For writing scalar steps at nodes
	xmlfile << "<Grid Name=\"mesh\" GridType=\"Uniform\">" << std::endl;
	Write_topology_data(Reader);
	Write_geometry_data(Reader);
	xmlfile << "<Time Value=\""+ std::to_string(time) + "\" />" << std::endl;
	xmlfile << "<Attribute Name=\""+data_name+"\" AttributeType=\"Scalar\" Center =\"Node\">" << std::endl;
	xmlfile << "<DataItem Dimensions=\"" + std::to_string(Reader.num_nodes) + " 1\" Format=\"XML\">" << std::endl;
	// Write data
	for (int i = 0; i < Reader.num_nodes; i++) {
		xmlfile << std::to_string(data[i]) << std::endl;
	}

	xmlfile << "</DataItem>" << std::endl;
	xmlfile << "</Attribute>" << std::endl;
	xmlfile << "</Grid>" << std::endl;
}
