#include <MishMesh/OBJ.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include <MishMesh/macros.h>

using namespace std;

vector<string> split(const string input, char delim = ' ') {
	stringstream ss(input);
	vector<string> parts;
	string item;
	while(getline(ss, item, delim)) {
		parts.push_back(item);
	}
	return parts;
}

/**
 * Write a mesh to OBJ using the precision given in the mesh traits.
 * @param filename The filename.
 * @param mesh The mesh
 * @return True when the mesh could be written.
 * @tparam MeshT An OpenMesh mesh type.
 * @note Vertex normals and halfedge texcoord2D (face vertex UVs) are supported.
 */
template<typename MeshT>
bool MishMesh::writeOBJ(const MeshT &mesh, const string filename, int precision) {
	ofstream ofs(filename);
	ofs << std::setprecision(precision) << std::fixed;
	if(!ofs.is_open()) {
		return false;
	}
	for(auto vh : mesh.vertices()) {
		auto v = mesh.point(vh);
		ofs << "v " << v[0] << " " << v[1] << " " << v[2] << endl;
	}
	if(mesh.has_face_normals()) {
		for(auto fh : mesh.faces()) {
			auto fn = mesh.normal(fh);
			ofs << "vn " << fn[0] << " " << fn[1] << " " << fn[2] << endl;
		}
	}
	if(mesh.has_vertex_texcoords2D()) {
		for(auto vh : mesh.vertices()) {
			auto uv = mesh.texcoord2D(vh);
			ofs << "vt " << uv[0] << " " << uv[1] << endl;
		}
	}

	for(auto fh : mesh.faces()) {
		ofs << "f";
		FOR_CFH(h_it, fh) {
			auto vh = mesh.to_vertex_handle(*h_it);
			ofs << " " << vh.idx() + 1;
			if(mesh.has_halfedge_texcoords2D()) {
				ofs << "/";
				ofs << vh.idx() + 1;
				if(mesh.has_face_normals()) {
					ofs << "/" << fh.idx() + 1;
				}
			} else if(mesh.has_face_normals()) {
				ofs << "//" << fh.idx() + 1;
			}
		}
		ofs << endl;
	}
	ofs.close();
	return true;
}

/**
 * Read a mesh to OBJ using the precision given in the mesh traits.
 * @param filename The filename.
 * @param mesh The mesh
 * @param Use the UV coordinates and per-face UV ids instead of the vertices of the mesh.
 *        This can be used to extract the UVs of the mesh.
 * @tparam MeshT An OpenMesh mesh type.
 * @return True when the mesh could be read.
 * @note Vertex normals and halfedge texcoord2D are currently not supported.
 */
template<typename MeshT>
bool MishMesh::readOBJ(MeshT &mesh, const string filename, bool build_uv_mesh) {
	ifstream ifs(filename);
	if(!ifs.is_open()) {
		return false;
	}

	vector<pair<double, double>> UVs;
	map<pair<size_t, short>, size_t> facevertex_uvs;

	mesh.clear();

	string line;
	while(!ifs.eof()) {
		getline(ifs, line);
		if(line.substr(0, 2) == "v ") {
			if(build_uv_mesh) continue; // We are building a mesh with the UVs as vertex coordinates
			vector<string> parts = split(line);
			array<double, 3> xyz;
			if(parts.size() == 4 || parts.size() == 7) {
				// v x y z [r g b]
				std::transform(parts.begin() + 1, parts.end(), xyz.begin(), [&](string part) { return std::stod(part); });
				mesh.add_vertex(typename MeshT::Point(xyz.data()));
			} else {
				cerr << "invalid vertex in OBJ file" << endl;
				return false;
			}
		} else if(line.substr(0, 3) == "vt ") {
			vector<string> parts = split(line);
			pair<double, double> uv;
			if(parts.size() == 3) {
				uv.first = std::stod(parts[1]);
				uv.second = std::stod(parts[2]);
			}
			UVs.push_back(uv);
			if(build_uv_mesh) {
				mesh.add_vertex({uv.first, uv.second, 0.0});
			}
		} else if(line[0] == 'f') {
			size_t t_idx = mesh.n_faces();
			std::vector<int> t;
			vector<string> idxs = split(line);
			// possible formats for a part:
			// vertex_id/uv_id
			// vertex_id/uv_id/normal_id
			// vertex_id//normal_id

			idxs.erase(idxs.begin()); // remove "f"
			for(auto idx : idxs) {
				auto parts = split(idx, '/');
				if(!build_uv_mesh) {
					t.push_back(std::stoi(parts[0]) - 1);
				} else {
					assert(parts.size() > 1);
					assert(parts[1] != "");
					t.push_back(std::stoi(parts[1]) - 1);
				}
			}
			std::vector<typename MeshT::VertexHandle> vhs(t.size());
			std::transform(t.begin(), t.end(), vhs.begin(), [&](int i){ return mesh.vertex_handle(i);});
			mesh.add_face(vhs.data(), vhs.size());
		}
	}
	return true;
}

template bool MishMesh::writeOBJ(const MishMesh::TriMesh &mesh, const string filename, int precision);
template bool MishMesh::writeOBJ(const MishMesh::PolyMesh &mesh, const string filename, int precision);
template bool MishMesh::readOBJ(MishMesh::TriMesh &mesh, const string filename, bool build_uv_mesh);
template bool MishMesh::readOBJ(MishMesh::PolyMesh &mesh, const string filename, bool build_uv_mesh);
