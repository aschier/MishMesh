#include <iostream>

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <MishMesh/utils.h>
#include <MishMesh/smoothing.h>
#include <MishMesh/split.h>
#include <MishMesh/dual_contouring.h>

#include <acc/bvh_tree.h>

#include "../../thirdparty/ProgramOptions.hxx"

typedef acc::BVHTree<unsigned int, OpenMesh::Vec3d> BVHTree;

void expand_bounding_box(MishMesh::BBox<> &bbox, double offset) {
   for(short j = 0; j < 3; j++) {
       bbox.ltf[j] -= offset * 1.8;
       bbox.rbn[j] += offset * 1.8;
   }
}

void adjust_resolution(const MishMesh::BBox<> &bbox, uint &resolutionX, uint &resolutionY, uint &resolutionZ, uint resolution, bool cubes) {
   if(cubes) {
       double min_length = bbox.min_side_length();
       double scaleX = (bbox.rbn[0] - bbox.ltf[0]) / min_length;
       double scaleY = (bbox.rbn[1] - bbox.ltf[1]) / min_length;
       double scaleZ = (bbox.rbn[2] - bbox.ltf[2]) / min_length;

       resolutionX = static_cast<int>(resolution * scaleX);
       resolutionY = static_cast<int>(resolution * scaleY);
       resolutionZ = static_cast<int>(resolution * scaleZ);
       if(resolutionX != resolution || resolutionY != resolution || resolutionZ != resolution) {
           std::cerr << "Automatic resolution for cubic cells: " << resolutionX << "x" << resolutionY << "x" << resolutionZ << std::endl;
       }
   }
}

MishMesh::PolyMesh get_largest_component(const MishMesh::PolyMesh &mesh) {
   auto components = MishMesh::split_connected_components(mesh);
   std::vector<double> diameters(components.size());
   std::transform(components.begin(), components.end(), diameters.begin(),
                  [&](auto &mesh) { return MishMesh::bounding_box(mesh).diameter(); });
   size_t idx = std::distance(diameters.begin(), std::max_element(diameters.begin(), diameters.end()));
   return components[idx];
}

std::shared_ptr<BVHTree> create_bvh_tree(const MishMesh::TriMesh& mesh) {
    std::vector<OpenMesh::Vec3d> vertices;
    vertices.reserve(mesh.n_vertices());
    for(const auto &vh : mesh.vertices()) {
        const auto &point = mesh.point(vh);
        vertices.push_back(point);
    }

    std::vector<unsigned int> faces;
    faces.reserve(mesh.n_faces() * 3);
    for(auto fh : mesh.faces()) {
        auto vhs = MishMesh::face_vertices(mesh, fh);
        faces.push_back(vhs[0].idx());
        faces.push_back(vhs[1].idx());
        faces.push_back(vhs[2].idx());
    }

    return BVHTree::create(faces, vertices);
}

struct MeshDistanceParam {
	std::shared_ptr<BVHTree> tree;
	double offset;
	const bool use_unsigned_distance = true;
};

/**
 * Calculate the offset distance and normal for a point relative to the mesh.
 *
 * @param p The point to calculate offset for
 * @param data Pointer to MeshDistanceParam containing BVH tree and offset parameters
 * @return Pair containing (offset distance, normal vector)
 */
std::pair<double, OpenMesh::Vec3d> distance_function(OpenMesh::Vec3d p, void *data) {
	MeshDistanceParam *param = static_cast<MeshDistanceParam *>(data);
	BVHTree::ClosestHit closest_hit;
	param->tree->closest_point(p, &closest_hit);
	auto vec = closest_hit.vertex - p;
	double sign = 1.0;
	if(!param->use_unsigned_distance && !closest_hit.front_face) {
		sign = -1.0;
	}
	return std::pair<double, OpenMesh::Vec3d>(sign * vec.norm() - param->offset, vec.normalized());
}

po::parser create_cmdline_parser() {
	po::parser parser;

	// clang-format off
	parser["help"]
		.abbreviation('h')
		.description("This help text")
		.callback([&]{
		std::cout << parser << '\n';
		exit(0);
	});
	parser["input"]
		.abbreviation('i')
		.type(po::string)
		.description("Input mesh file path");
	parser["output"]
		.abbreviation('o')
		.type(po::string)
		.description("Output mesh file path");
	parser["offset"]
		.abbreviation('O')
		.type(po::f32)
		.description("Offset distance for mesh extraction (default: 3.0)")
		.fallback(3.0);
	parser["absolute"]
		.abbreviation('a')
		.description("Use absolute offset value instead of percentage of bounding box diagonal");
	parser["resolution"]
		.abbreviation('r')
		.type(po::i32)
		.description("Number of voxels to sample in each dimension (default: 60)")
		.fallback(60);
	parser["cubes"]
		.abbreviation('c')
		.description("Adjust resolution for cubic grid cells");
	parser["smoothing"]
		.abbreviation('s')
		.type(po::i32)
		.description("Number of smoothing steps to apply (default: 0)")
		.fallback(0);
	parser["triangulate"]
		.abbreviation('t')
		.description("Triangulate the output mesh if it's not already triangulated");
	parser["unsigned-distance"]
		.abbreviation('u')
		.description("Use unsigned distance for two-sided offset mesh (ignores front/back face orientation)");
	parser["largest-component"]
		.abbreviation('L')
		.description("Keep only the largest connected component in the output");
	// clang-format on
	return parser;
}

int main(int argc, char **argv) {
	po::parser parser = create_cmdline_parser();
	parser.parse(argc, argv);

	if(argc == 1) {
		std::cout << parser << '\n';
		exit(1);
	}
	if(!parser["input"].available()) {
		std::cerr << "Error: No input file was provided. Please specify an input mesh using the -i or --input option." << std::endl;
		exit(1);
	}
	if(!parser["output"].available()) {
		std::cerr << "Error: No output filename was provided. Please specify an output mesh using the -o or --output option." << std::endl;
		exit(1);
	}

	std::string filename(parser["input"].get().string);
	std::string output_filename(parser["output"].get().string);

	double offset = parser["offset"].get().f32;
	uint resolution = parser["resolution"].get().i32;
	bool offset_is_absolute = parser["absolute"].available();
	uint smoothing_steps = parser["smoothing"].get().i32;
	bool triangulate = parser["triangulate"].available();
	bool cubes = parser["cubes"].available();
	bool use_unsigned_distance = parser["unsigned-distance"].available();
	bool only_largest_component = parser["largest-component"].available();

	MishMesh::TriMesh input_mesh;
	bool success = OpenMesh::IO::read_mesh(input_mesh, filename);
	if(!success) {
		std::cerr << "Error: Failed to read the input mesh file '" << filename << "'. Please check that the file exists and is a valid mesh." << std::endl;
		exit(1);
	}

	// Use the same resolution for all dimensions initially
	uint resolutionX = resolution;
	uint resolutionY = resolution;
	uint resolutionZ = resolution;

	MishMesh::BBox<> bbox = MishMesh::bounding_box(input_mesh);

	adjust_resolution(bbox, resolutionX, resolutionY, resolutionZ, resolution, cubes);

	if(!offset_is_absolute) {
		offset = offset / 100 * (bbox.rbn - bbox.ltf).norm();
	}

	expand_bounding_box(bbox, offset);

	auto bvh_tree = create_bvh_tree(input_mesh);
	MeshDistanceParam param{bvh_tree, offset, use_unsigned_distance};
	MishMesh::PolyMesh mesh = MishMesh::Meshing::dual_contouring(distance_function, &param, bbox, {resolutionX, resolutionY, resolutionZ},
	                                                             MishMesh::Meshing::VertexFit::Average);

	if(only_largest_component) {
		mesh = get_largest_component(mesh);
	}

	MishMesh::smooth_mesh(mesh, smoothing_steps);
	if(triangulate) {
		mesh.triangulate();
	}

	OpenMesh::IO::write_mesh(mesh, output_filename);

	return 0;
}
