//
#include "weno.ispc"
//
export void
#ifdef USE_FLOAT
	ispc_interpolate_trilinear_float3(
#elif USE_DOUBLE
	ispc_interpolate_trilinear_double3(
#endif
	uniform uint8 head0, uniform uint32 size,
	uniform const T value_cache[],
	uniform const T x_array[],
	uniform const T y_array[],
	uniform const T z_array[],
	uniform T result[] ) {
	//
	uniform const uint8 w = 2;
	uniform const uint8 ww = w*w;
	uniform const uint8 head1 = (head0+1)%w;
	foreach (n = 0 ... size) {
		//
		const T x = x_array[n];
		const T y = y_array[n];
		const T z = z_array[n];
		//
		const T v000 = value_cache[0+w*0+head0*ww];
		const T v100 = value_cache[1+w*0+head0*ww];
		const T v010 = value_cache[0+w*1+head0*ww];
		const T v110 = value_cache[1+w*1+head0*ww];
		const T v001 = value_cache[0+w*0+head1*ww];
		const T v101 = value_cache[1+w*0+head1*ww];
		const T v011 = value_cache[0+w*1+head1*ww];
		const T v111 = value_cache[1+w*1+head1*ww];
		//
		result[n] = 
			(1.0-x)*(1.0-y)*(1-z)*v000 + (x)*(1.0-y)*(1-z)*v100 +
			(1.0-x)*(y)*(1-z)*v010 + (x)*(y)*(1-z)*v110 +
			(1.0-x)*(1.0-y)*(z)*v001 + (x)*(1.0-y)*(z)*v101 +
			(1.0-x)*(y)*(z)*v011 + (x)*(y)*(z)*v111;
	}
}
//
export void
#ifdef USE_FLOAT
	ispc_interpolate_weno4_float3(
#elif USE_DOUBLE
	ispc_interpolate_weno4_double3(
#endif
	uniform uint8 head0, uniform uint32 size,
	uniform const T value_cache[],
	uniform const T x_array[],
	uniform const T y_array[],
	uniform const T z_array[],
	uniform T result[] ) {
	//
	uniform const uint8 w = 4;
	uniform const uint8 ww = w*w;
	foreach (n = 0 ... size) {
		//
		const T x = x_array[n];
		const T y = y_array[n];
		const T z = z_array[n];
		//
		T tmp_k[4], tmp_x[4];
		for( uniform int k=0; k<4; ++k ) {
			uniform const uint8 head = (head0+k)%w;
			uniform T S[2];
			weno_interp4_uniform_precompute(S,value_cache+0*w+head*ww);
			tmp_x[0] = weno_interp4_uniform(x,S,value_cache+0*w+head*ww);
			tmp_x[1] = weno_interp4_uniform(x,S,value_cache+1*w+head*ww);
			tmp_x[2] = weno_interp4_uniform(x,S,value_cache+2*w+head*ww);
			tmp_x[3] = weno_interp4_uniform(x,S,value_cache+3*w+head*ww);
			tmp_k[k] = weno_interp4(y,tmp_x);
		}
		result[n] = weno_interp4(z,tmp_k);
	}
}
//
export void
#ifdef USE_FLOAT
	ispc_interpolate_weno6_float3(
#elif USE_DOUBLE
	ispc_interpolate_weno6_double3(
#endif
	uniform uint8 head0, uniform uint32 size,
	uniform const T value_cache[],
	uniform const T x_array[],
	uniform const T y_array[],
	uniform const T z_array[],
	uniform T result[] ) {
	//
	uniform const uint8 w = 6;
	uniform const uint8 ww = w*w;
	foreach (n = 0 ... size) {
		//
		const T x = x_array[n];
		const T y = y_array[n];
		const T z = z_array[n];
		//
		T tmp_k[6], tmp_x[6];
		for( uniform int k=0; k<6; ++k ) {
			uniform const uint8 head = (head0+k)%w;
			uniform T S[3];
			weno_interp6_uniform_precompute(S,value_cache+0*w+head*ww);
			tmp_x[0] = weno_interp6_uniform(x,S,value_cache+0*w+head*ww);
			tmp_x[1] = weno_interp6_uniform(x,S,value_cache+1*w+head*ww);
			tmp_x[2] = weno_interp6_uniform(x,S,value_cache+2*w+head*ww);
			tmp_x[3] = weno_interp6_uniform(x,S,value_cache+3*w+head*ww);
			tmp_x[4] = weno_interp6_uniform(x,S,value_cache+4*w+head*ww);
			tmp_x[5] = weno_interp6_uniform(x,S,value_cache+5*w+head*ww);
			tmp_k[k] = weno_interp6(y,tmp_x);
		}
		result[n] = weno_interp6(z,tmp_k);
	}
}