#include <MishMesh/VTK.h>

#include <fstream>
#include <iomanip>

#include <MishMesh/utils.h>

using TriMesh = MishMesh::TriMesh;

void MishMesh::writeVTK(const TriMesh &mesh, const std::string filename,
              VertexValueProperty prop_vertex_value,
              VertexVectorProperty prop_vertex_vector,
              FaceValueProperty prop_face_value,
              FaceVectorProperty prop_face_vector) {

	std::ofstream ofs(filename);

	// Write header
	ofs << std::fixed << std::setprecision(6);
	ofs << "# vtk DataFile Version 2.0" << std::endl;
	ofs << "Mesh" << std::endl;
	ofs << "ASCII" << std::endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	ofs << "POINTS " << mesh.n_vertices() << " float" << std::endl;

	// Write vertices
	for(auto vh : mesh.vertices()) {
		auto point = mesh.point(vh);
		ofs << point[0] << " " << point[1] << " " << point[2] << " " << std::endl;
	}

	// Write faces
	ofs << "CELLS " << mesh.n_faces() << " " << mesh.n_faces() * (3 + 1) << std::endl;
	for(auto fh : mesh.faces()) {
		auto vhs = MishMesh::face_vertices(mesh, fh);
		ofs << "3 " << vhs[0].idx() << " " << vhs[1].idx() << " " << vhs[2].idx() << " " << std::endl;
	}
	// Triangle cells
	ofs << "CELL_TYPES " << mesh.n_faces() << std::endl;
	for(uint i = 0; i < mesh.n_faces(); i++) {
		ofs << 5 << std::endl;
	}

	// Write the vertex values and vertex vectors
	if(prop_vertex_value.is_valid() || prop_vertex_vector.is_valid()) {
		ofs << "POINT_DATA " << mesh.n_vertices() << std::endl;
	}
	if(prop_vertex_value.is_valid()) {
		ofs << "SCALARS "
		    << "point_values"
		    << " float" << std::endl;
		ofs << "LOOKUP_TABLE default" << std::endl;
		for(auto vh : mesh.vertices()) {
			ofs << mesh.property(prop_vertex_value, vh) << std::endl;
		}
	}
	if(prop_vertex_vector.is_valid()) {
		ofs << "VECTORS "
		    << "point_vector"
		    << " float" << std::endl;
		for(auto vh : mesh.vertices()) {
			auto v = mesh.property(prop_vertex_vector, vh);
			ofs << v[0] << " " << v[1] << " " << v[2] << std::endl;
		}
	}

	// Write the face values and face vectors
	if(prop_face_value.is_valid() || prop_face_vector.is_valid()) {
		ofs << "CELL_DATA " << mesh.n_faces() << std::endl;
	}
	if(prop_face_value.is_valid()) {
		ofs << "SCALARS "
		    << "face_value"
		    << " float" << std::endl;
		ofs << "LOOKUP_TABLE default" << std::endl;
		for(auto fh : mesh.faces()) {
			ofs << mesh.property(prop_face_value, fh) << std::endl;
		}
	}
	if(prop_face_vector.is_valid()) {
		ofs << "VECTORS "
		    << "face_vectors"
		    << " float" << std::endl;
		for(auto fh : mesh.faces()) {
			auto v = mesh.property(prop_face_vector, fh);
			ofs << v[0] << " " << v[1] << " " << v[2] << std::endl;
		}
	}
	ofs.close();
}
