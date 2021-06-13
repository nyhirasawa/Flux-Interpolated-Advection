#include "vdbexport.h"

#include <openvdb/openvdb.h>
#include <cmath>

#include "define_float_type.h"
#include "utils.h"
//void export2vdb(
//	const char *path,
//	unsigned width, unsigned height, unsigned depth,
//	MY_FLOAT_TYPE ox, MY_FLOAT_TYPE oy, MY_FLOAT_TYPE oz, MY_FLOAT_TYPE dx,
//	const MY_FLOAT_TYPE *data
//) {
	//
//	openvdb::initialize();
//    openvdb::DoubleGrid::Ptr grid = openvdb::DoubleGrid::create(/*background value=*/0.0);
	//
//	openvdb::math::Transform::Ptr transform =
//		openvdb::math::Transform::createLinearTransform(/*voxel size=*/dx);
//	transform->postTranslate(openvdb::math::Vec3d(ox,oy,oz));
//	grid->setTransform(transform);
	//
//	openvdb::DoubleGrid::Accessor accessor = grid->getAccessor();
//	const size_t xy_count = width*height;
//	for(size_t n=0; n<xy_count*depth; n++) {
//		const unsigned k = n / xy_count;
//		const size_t m (n % xy_count);
//		const unsigned i (m % width);
//		const unsigned j (m / width);
//		const openvdb::Coord xyz(i,j,k);
//		const MY_FLOAT_TYPE value = data[i+j*width+k*xy_count];
//		if( value ) {
//			accessor.setValue(xyz,value);
//		}
//	}
//	openvdb::io::File file(path);
//	openvdb::GridPtrVec grids;
//	grids.push_back(grid);
//	file.write(grids);
//	file.close();
//}

void export2vdb(
	const char *path,
	unsigned width, unsigned height, unsigned depth,
	MY_FLOAT_TYPE ox, MY_FLOAT_TYPE oy, MY_FLOAT_TYPE oz, MY_FLOAT_TYPE dx,
	const MY_FLOAT_TYPE *data
) {
	//
	openvdb::initialize();
	openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0);
	//
	openvdb::math::Transform::Ptr transform =
		openvdb::math::Transform::createLinearTransform(dx);
	transform->postTranslate(openvdb::math::Vec3d(ox,oy,oz));
	grid->setName("density");
	grid->setTransform(transform);
	//
	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
	for(int ix=0;ix<width; ++ix){
		for(int iy=0;iy<height; ++iy){
			for(int iz=0;iz<depth; ++iz){
				const openvdb::Coord xyz(ix,iy,iz);
//				const float value = data[smoke_simulation::get_voxel_center_index_3D(width-1-ix, iy, iz, width, height, depth)];
				const float value = data[smoke_simulation::get_voxel_center_index_3D(ix, iy, iz, width, height, depth)];

				accessor.setValue(xyz, value);
/*
				if( value < 0.0) {
					accessor.setValue(xyz, -value);
				}
				else{
					accessor.setValue(xyz, 0.0);
				}
*/
			}
		}
	}
/*
	const size_t xy_count = width*height;
	for(size_t n=0; n<xy_count*depth; n++) {
		const unsigned k = n / xy_count;
		const size_t m (n % xy_count);
		const unsigned i (m % width);
		const unsigned j (m / width);
		const openvdb::Coord xyz(i,j,k);
		const float value = data[k+j*width+i*xy_count];
		if( value ) {
			accessor.setValue(xyz,value);
		}
	}
*/
	openvdb::io::File file(path);
	openvdb::GridPtrVec grids;
	grids.push_back(grid);
	file.write(grids);
	file.close();
}
