#include <OpenMesh/Core/IO/MeshIO.hh>

#include <MishMesh/PolyMesh.h>
#include <MishMesh/dual_contouring.h>

struct SphereParameter {
	OpenMesh::Vec3d center;
	double radius;
};

std::pair<double, OpenMesh::Vec3d> sphere(const OpenMesh::Vec3d p, void *data = nullptr) {
	SphereParameter *sphere_param = (SphereParameter *)data;
	double distance = (p - sphere_param->center).sqrnorm() - (sphere_param->radius * sphere_param->radius);
	auto normal = (p - sphere_param->center).normalized();
	if(normal.norm() == 0) {
		// When a grid point lies on the center of the sphere, there is no useful normal.
		normal = {1.0, 0, 0};
	}
	assert(!std::isnan(normal.norm()));
	return std::make_pair(distance, normal);
}

struct CSGSphereParameter {
	SphereParameter sphere_param1, sphere_param2;
};

struct PlaneParameter {
	short dim;
	double coordinate;
};

std::pair<double, OpenMesh::Vec3d> plane(OpenMesh::Vec3d p, void *data) {
	PlaneParameter *param = static_cast<PlaneParameter *>(data);
	return {
	    param->coordinate - p[param->dim],
	    OpenMesh::Vec3d{param->dim == 0 ? 1 : 0, param->dim == 1 ? 1 : 0, param->dim == 2 ? 1 : 0}};
}

enum CSGMode {
	MIN,
	MAX
};

struct CSGParameter {
	std::pair<double, OpenMesh::Vec3d> (*distance_function1)(OpenMesh::Vec3d, void *);
	std::pair<double, OpenMesh::Vec3d> (*distance_function2)(OpenMesh::Vec3d, void *);
	void *param1;
	void *param2;
	CSGMode mode = MIN;
};

std::pair<double, OpenMesh::Vec3d> csg(OpenMesh::Vec3d p, void *data) {
	CSGParameter *csg_param = static_cast<CSGParameter *>(data);
	auto result1 = csg_param->distance_function1(p, csg_param->param1);
	auto result2 = csg_param->distance_function2(p, csg_param->param2);
	if((result1.first < result2.first && csg_param->mode == MIN) ||
	   (result1.first >= result2.first && csg_param->mode == MAX)) {
		return result1;
	} else {
		return result2;
	}
}

int main(int argc, char **argv) {
	SphereParameter sphere_param1{{0, 0, 0}, 0.9};
	SphereParameter sphere_param2{{0, 1.5, 0}, 0.9};

	CSGParameter param{
	    sphere,
	    sphere,
	    &sphere_param1, &sphere_param2,
	    CSGMode::MIN};

	std::cerr << "Generating merging spheres with averaging vertex fitting (good for round shapes): spheres_average.obj" << std::endl;
	MishMesh::PolyMesh mesh = MishMesh::Meshing::dual_contouring(
	    csg, &param,
	    MishMesh::BBox<OpenMesh::Vec3d>(
	        OpenMesh::Vec3d{-1, -1, -1},
	        OpenMesh::Vec3d{1, 2.5, 1}),
	    OpenMesh::Vec3i{50, 50, 50},
	    MishMesh::Meshing::VertexFit::Average);
	OpenMesh::IO::write_mesh(mesh, "spheres_average.obj");

#ifdef HAS_EIGEN
	std::cerr << "Generating merging spheres with SVD vertex fitting (with vertex fitting artifacts): spheres_svd.obj" << std::endl;
	mesh = MishMesh::Meshing::dual_contouring(
	    csg, &param,
	    MishMesh::BBox<OpenMesh::Vec3d>(
	        OpenMesh::Vec3d{-1, -1, -1},
	        OpenMesh::Vec3d{1, 2.5, 1}),
	    OpenMesh::Vec3i{50, 50, 50},
	    MishMesh::Meshing::VertexFit::SVD);
	OpenMesh::IO::write_mesh(mesh, "spheres_svd.obj");
#endif

	PlaneParameter plane_param1{0, 0.1};
	PlaneParameter plane_param2{1, 0.1};
	param = {
	    plane, plane,
	    &plane_param1, &plane_param2,
	    CSGMode::MAX};

	std::cerr << "Generating planes with averaging vertex fitting (losing sharp features): planes_average.obj" << std::endl;
	mesh = MishMesh::Meshing::dual_contouring(
	    csg, &param,
	    MishMesh::BBox<OpenMesh::Vec3d>(
	        OpenMesh::Vec3d{0, 0, 0},
	        OpenMesh::Vec3d{1, 1, 1}),
	    OpenMesh::Vec3i{50, 50, 50},
	    MishMesh::Meshing::VertexFit::Average);
	OpenMesh::IO::write_mesh(mesh, "planes_average.obj");

#ifdef HAS_EIGEN
	std::cerr << "Generating planes with SVD vertex fitting (preserving sharp features): planes_svd.obj" << std::endl;
	mesh = MishMesh::Meshing::dual_contouring(
	    csg, &param,
	    MishMesh::BBox<OpenMesh::Vec3d>(
	        OpenMesh::Vec3d{0, 0, 0},
	        OpenMesh::Vec3d{1, 1, 1}),
	    OpenMesh::Vec3i{50, 50, 50},
	    MishMesh::Meshing::VertexFit::SVD);
	OpenMesh::IO::write_mesh(mesh, "planes_svd.obj");
#endif
}
