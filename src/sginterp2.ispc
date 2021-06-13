//
#include "weno.ispc"
//
export void
#ifdef USE_FLOAT
	ispc_interpolate_bilinear_float2(
#elif USE_DOUBLE
	ispc_interpolate_bilinear_double2(
#endif
	uniform uint8 head0, uniform uint32 size,
	uniform const T value_cache[],
	uniform const T x_array[],
	uniform const T y_array[],
	uniform T result[] ) {
	//
	uniform const uint8 w = 2;
	uniform const uint8 head1 = (head0+1)%w;
	foreach (n = 0 ... size) {
		//
		const T x = x_array[n];
		const T y = y_array[n];
		//
		const T v00 = value_cache[w*head0+0];
		const T v10 = value_cache[w*head0+1];
		const T v01 = value_cache[w*head1+0];
		const T v11 = value_cache[w*head1+1];
		//
		result[n] = 
			(1.0-x)*(1.0-y)*v00 + (x)*(1.0-y)*v10 +
			(1.0-x)*(y)*v01 + (x)*(y)*v11;
	}
}
//
export void
#ifdef USE_FLOAT
	ispc_interpolate_weno4_float2(
#elif USE_DOUBLE
	ispc_interpolate_weno4_double2(
#endif
	uniform uint8 head0, uniform uint32 size,
	uniform const T value_cache[],
	uniform const T x_array[],
	uniform const T y_array[],
	uniform T result[] ) {
	//
	uniform const uint8 w = 4;
	foreach (n = 0 ... size) {
		//
		const T x = x_array[n];
		const T y = y_array[n];
		//
		T tmp_x[4];
		for( uniform uint8 m=0; m<4; ++m ) {
			const uniform uint8 head = (head0+m)%w;
			uniform T S[2];
			weno_interp4_uniform_precompute(S,value_cache+head*w);
			tmp_x[m] = weno_interp4_uniform(x,S,value_cache+head*w);
		}
		result[n] = weno_interp4(y,tmp_x);
	}
}
//
export void
#ifdef USE_FLOAT
	ispc_interpolate_weno6_float2(
#elif USE_DOUBLE
	ispc_interpolate_weno6_double2(
#endif
	uniform uint8 head0, uniform uint32 size,
	uniform const T value_cache[],
	uniform const T x_array[],
	uniform const T y_array[],
	uniform T result[] ) {
	//
	uniform const uint8 w = 6;
	foreach (n = 0 ... size) {
		//
		const T x = x_array[n];
		const T y = y_array[n];
		//
		T tmp_x[6];
		for( uniform uint8 m=0; m<6; ++m ) {
			const uniform uint8 head = (head0+m)%w;
			uniform T S[3];
			weno_interp6_uniform_precompute(S,value_cache+head*w);
			tmp_x[m] = weno_interp6_uniform(x,S,value_cache+head*w);
		}
		result[n] = weno_interp6(y,tmp_x);
	}
}