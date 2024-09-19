#pragma once
#include <iostream>
#include <fstream>
#include "Mesh.h"

class XML_Writer
{
private:
	std::ofstream xmlfile;
public:
	XML_Writer(std::string path_storage, std::string filename);
	~XML_Writer();
	void Write_topology_data(Mesh& Reader);
	void Write_geometry_data(Mesh& Reader);
	void Write_time_step(Mesh& Reader, double time, std::string data_name, std::vector<double>& data);
};

