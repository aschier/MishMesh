#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>

#include <MishMesh/TriMesh.h>
#include <MishMesh/sampling.h>
#include <MishMesh/poisson_disk_sampling.h>

#include <array>
#include <chrono>

int main(int argc, char **argv) {
	typedef std::chrono::high_resolution_clock clock;
	typedef std::chrono::duration<double, std::milli> duration;

	MishMesh::TriMesh mesh;
	std::array<OpenMesh::Vec3d, 3> triangle_points{
		OpenMesh::Vec3d{0.0, 0.0, 0.0},
		OpenMesh::Vec3d{1.0, 0.0, 0.0},
		OpenMesh::Vec3d{0.0, 1.0, 0.0}
	};

	std::cerr << "Sampling uniform inside a triangle:" << std::endl;
	auto timepoint = clock::now();
	for(int i = 0; i < 500; i++) {
		auto p = MishMesh::uniform_random_triangle_point(triangle_points);
		mesh.add_vertex(p);
	}
	duration time = clock::now() - timepoint;
	std::cerr << "Sampled " << mesh.n_vertices() << " points in " << time.count() << "ms." << std::endl << std::endl;
	OpenMesh::IO::write_mesh(mesh, "uniform_triangle.obj");

	// Uniform annulus sampling
	mesh.clear();
	std::cerr << "Sampling uniform inside an annulus:" << std::endl;
	timepoint = clock::now();
	for(int i = 0; i < 500; i++) {
		auto p = MishMesh::uniform_random_annulus_point({0, 0}, 1.0, 0.5);
		mesh.add_vertex({p[0], p[1], 0.0});
	}
	time = clock::now() - timepoint;
	std::cerr << "Sampled " << mesh.n_vertices() << " points in " << time.count() << "ms." << std::endl << std::endl;
	OpenMesh::IO::write_mesh(mesh, "uniform_annulus.obj");

	// Poisson sampling
	mesh.clear();
	std::cerr << "Poisson sampling inside a triangle:" << std::endl;
	for(int i = 1; i < 15; i++) {
		std::array<OpenMesh::Vec2d, 3> triangle_points2D{
			OpenMesh::Vec2d{0.1 + 1.1*i, 0.0},
			OpenMesh::Vec2d{1.0 + 1.1*i, 0.2},
			OpenMesh::Vec2d{0.7 + 1.1*i, 0.9}
		};
		double min_dist = std::pow(10, -i * 0.15);
		timepoint = clock::now();
		auto points = MishMesh::poisson_disk_sampling(triangle_points2D, min_dist, 0.5 * min_dist / sqrt(2));
		time = clock::now() - timepoint;
		std::cerr << "Sampled " << points.size() << " points in " << time.count() << "ms." << std::endl;
		for(auto p : points){
			mesh.add_vertex({p[0], p[1], 0.0});
		}
	}
	std::cerr << std::endl;
	OpenMesh::IO::write_mesh(mesh, "poisson_sampling.obj");
}
