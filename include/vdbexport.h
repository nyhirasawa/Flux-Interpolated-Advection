#pragma once

#include <openvdb/openvdb.h>
#include <cmath>

#include "define_float_type.h"

//namespace smoke_simulation{

//
/*
void export2vdb(
	const char *path,
	unsigned width, unsigned height, unsigned depth,
	MY_FLOAT_TYPE ox, MY_FLOAT_TYPE oy, MY_FLOAT_TYPE oz, MY_FLOAT_TYPE dx,
	const MY_FLOAT_TYPE *data
);
*/

void export2vdb(
	const char *path,
	unsigned width, unsigned height, unsigned depth,
	MY_FLOAT_TYPE ox, MY_FLOAT_TYPE oy, MY_FLOAT_TYPE oz, MY_FLOAT_TYPE dx,
	const MY_FLOAT_TYPE *data
);

//}
