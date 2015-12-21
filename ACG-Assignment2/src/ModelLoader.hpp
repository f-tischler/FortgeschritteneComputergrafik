#ifndef ModelLoader_H_included
#define ModelLoader_H_included

#include <memory>
#include <iostream>

#include "TriangleMesh.hpp"
#include "TinyObjLoader.hpp"

inline std::vector<std::shared_ptr<TriangleMesh>> loadObj(const std::string& filename, Color emmision, Color color, Refl_t reflType)
{
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err; 

	tinyobj::LoadObj(shapes, materials, err, filename.c_str());

	if (!err.empty()) {
		std::cerr << err << std::endl;
		exit(1);
	}

	std::vector<std::shared_ptr<TriangleMesh>> result;

	for (auto it = shapes.begin(); it != shapes.end(); ++it)
	{
		auto trimesh = std::make_shared<TriangleMesh>(emmision, color, reflType);

		const auto& shape = *it;
		for (size_t i = 0; i < shape.mesh.indices.size(); i++)
		{
			auto x = shape.mesh.positions[shape.mesh.indices[i]*3];
			auto y = shape.mesh.positions[shape.mesh.indices[i]*3+1];
			auto z = shape.mesh.positions[shape.mesh.indices[i]*3+2];

			trimesh->addVertex(Vertex(x,y,z));
		}


		result.push_back(std::move(trimesh));
	}

	return result;
}

#endif