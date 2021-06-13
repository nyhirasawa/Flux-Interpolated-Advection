//
#include "sginterp3.h"
#include "weno.h"
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <algorithm>
//
extern "C" {
	//
	void ispc_interpolate_trilinear_float3(
		unsigned char head, unsigned size,
		const float value_cache[],
		const float x_array[],
		const float y_array[],
		const float z_array[],
		float result[] );
	//
	void ispc_interpolate_trilinear_double3(
		unsigned char head, unsigned size,
		const double value_cache[],
		const double x_array[],
		const double y_array[],
		const double z_array[],
		double result[] );
	//
	void ispc_interpolate_weno4_float3(
		unsigned char head, unsigned size,
		const float value_cache[],
		const float x_array[],
		const float y_array[],
		const float z_array[],
		float result[] );
	//
	void ispc_interpolate_weno4_double3(
		unsigned char head, unsigned size,
		const double value_cache[],
		const double x_array[],
		const double y_array[],
		const double z_array[],
		double result[] );
	//
	void ispc_interpolate_weno6_float3(
		unsigned char head, unsigned size,
		const float value_cache[],
		const float x_array[],
		const float y_array[],
		const float z_array[],
		float result[] );
	//
	void ispc_interpolate_weno6_double3(
		unsigned char head, unsigned size,
		const double value_cache[],
		const double x_array[],
		const double y_array[],
		const double z_array[],
		double result[] );
}
//
template <class T>
static T clamp( T v, T min, T max ) {
	return std::max(min,std::min(v,max));
}
//
template <class T>
static void interpolate_trilinear(
	unsigned char head0, unsigned size, bool ispc,
	const T *value_cache,
	const T *x_array,
	const T *y_array,
	const T *z_array,
	T *result ) {
	//
	if( ispc ) {
		if( std::is_same<T,float>::value ) {
			ispc_interpolate_trilinear_float3(head0,size,
				(const float *)value_cache,
				(const float *)x_array,(const float *)y_array,(const float *)z_array,
				(float *)result);
		} else if( std::is_same<T,double>::value ) {
			ispc_interpolate_trilinear_double3(head0,size,
				(const double *)value_cache,
				(const double *)x_array,(const double *)y_array,(const double *)z_array,
				(double *)result);
		}
	} else {
		//
		const unsigned char w (2);
		const unsigned char ww (4);
		const unsigned char head1 ((head0+1)%w);
		for( unsigned n=0; n<size; ++n ) {
			//
			const T x (x_array[n]), y (y_array[n]), z (z_array[n]);
			//
			const T &v000 = value_cache[0+w*0+head0*ww];
			const T &v100 = value_cache[1+w*0+head0*ww];
			const T &v010 = value_cache[0+w*1+head0*ww];
			const T &v110 = value_cache[1+w*1+head0*ww];
			const T &v001 = value_cache[0+w*0+head1*ww];
			const T &v101 = value_cache[1+w*0+head1*ww];
			const T &v011 = value_cache[0+w*1+head1*ww];
			const T &v111 = value_cache[1+w*1+head1*ww];
			//
			result[n] = 
				(1.0-x)*(1.0-y)*(1-z)*v000 + (x)*(1.0-y)*(1-z)*v100 +
				(1.0-x)*(y)*(1-z)*v010 + (x)*(y)*(1-z)*v110 +
				(1.0-x)*(1.0-y)*(z)*v001 + (x)*(1.0-y)*(z)*v101 +
				(1.0-x)*(y)*(z)*v011 + (x)*(y)*(z)*v111;
		}
	}
}
//
template <class T>
static void interpolate_WENO4(
	unsigned char head0, unsigned size, bool ispc,
	const T *value_cache,
	const T *x_array,
	const T *y_array,
	const T *z_array,
	T *result ) {
	//
	if( ispc ) {
		if( std::is_same<T,float>::value ) {
			ispc_interpolate_weno4_float3(head0,size,
				(const float *)value_cache,
				(const float *)x_array,(const float *)y_array,(const float *)z_array,
				(float *)result);
		} else if( std::is_same<T,double>::value ) {
			ispc_interpolate_weno4_double3(head0,size,
				(const double *)value_cache,
				(const double *)x_array,(const double *)y_array,(const double *)z_array,
				(double *)result);
		}
	} else {
		//
		const unsigned char w (4);
		const unsigned char ww (16);
		for( unsigned n=0; n<size; ++n ) {
			//
			const T x (x_array[n]), y (y_array[n]), z (z_array[n]);
			//
			T tmp_k[4];
			for( int k=0; k<4; ++k ) {
				T tmp_x[4];
				const unsigned char head = (head0+k)%w;
				tmp_x[0] = WENO<T>::interp4(x,value_cache+0*w+head*ww);
				tmp_x[1] = WENO<T>::interp4(x,value_cache+1*w+head*ww);
				tmp_x[2] = WENO<T>::interp4(x,value_cache+2*w+head*ww);
				tmp_x[3] = WENO<T>::interp4(x,value_cache+3*w+head*ww);
				tmp_k[k] = WENO<T>::interp4(y,tmp_x);
			}
			result[n] = WENO<T>::interp4(z,tmp_k);
		}
	}
}
//
template <class T>
static void interpolate_WENO6(
	unsigned char head0, unsigned size, bool ispc,
	const T *value_cache,
	const T *x_array,
	const T *y_array,
	const T *z_array,
	T *result ) {
	//
	if( ispc ) {
		if( std::is_same<T,float>::value ) {
			ispc_interpolate_weno6_float3(head0,size,
				(const float *)value_cache,
				(const float *)x_array,(const float *)y_array,(const float *)z_array,
				(float *)result);
		} else if( std::is_same<T,double>::value ) {
			ispc_interpolate_weno6_double3(head0,size,
				(const double *)value_cache,
				(const double *)x_array,(const double *)y_array,(const double *)z_array,
				(double *)result);
		}
	} else {
		//
		const unsigned char w (6);
		const unsigned char ww (36);
		for( unsigned n=0; n<size; ++n ) {
			//
			const T x (x_array[n]), y (y_array[n]), z (z_array[n]);
			//
			T tmp_k[6];
			for( int k=0; k<6; ++k ) {
				T tmp_x[6];
				const unsigned char head = (head0+k)%w;
				tmp_x[0] = WENO<T>::interp6(x,value_cache+0*w+head*ww);
				tmp_x[1] = WENO<T>::interp6(x,value_cache+1*w+head*ww);
				tmp_x[2] = WENO<T>::interp6(x,value_cache+2*w+head*ww);
				tmp_x[3] = WENO<T>::interp6(x,value_cache+3*w+head*ww);
				tmp_x[4] = WENO<T>::interp6(x,value_cache+4*w+head*ww);
				tmp_x[5] = WENO<T>::interp6(x,value_cache+5*w+head*ww);
				tmp_k[k] = WENO<T>::interp6(y,tmp_x);
			}
			result[n] = WENO<T>::interp6(z,tmp_k);
		}
	}
}
//
template <class T>
std::vector<T> sginterp3::interpolate(
		const std::vector<point3<T> > &positions,
		unsigned nx, unsigned ny, unsigned nz, T dx,
		const T *volume,
		filter scheme,
		int flag
) {
	//
	std::vector<T> results(positions.size());
	//
	unsigned char w (0), ww (0);
	const unsigned nxy = nx * ny;
	const unsigned cell_xy = (nx-1)*(ny-1);
	const bool ispc = flag & backend::simd;
	//
	switch( scheme ) {
		case filter::trilinear:
			w = 2;
			break;
		case filter::WENO4:
			w = 4;
			break;
		case filter::WENO6:
			w = 6;
			break;
		default:
			printf("Undefined scheme! (%d)\n",(int)scheme);
			exit(1);
	}
	ww = w*w;
	//
	auto load_initial_values = [&]( int i, int j, int k, int w, T *value_cache ) {
		for( int ii=0; ii<w; ++ii ) {
			for( int jj=0; jj<w; ++jj ) {
				for( int kk=0; kk<w; ++kk ) {
					const int fi = clamp(i+ii-w/2+1,0,(int)nx-1);
					const int fj = clamp(j+jj-w/2+1,0,(int)ny-1);
					const int fk = clamp(k+kk-w/2+1,0,(int)nz-1);
					value_cache[ii+jj*w+kk*ww] = volume[fk*nxy+fj*nx+fi];
				}
			}
		}
	};
	//
	if( flag & backend::presort ) {
		//
		std::vector<std::vector<size_t> > hash((nx-1)*(ny-1)*(nz-1));
		std::vector<std::mutex> mutex((nx-1)*(ny-1)*(nz-1));
		//
		printf( "Building hash table...\n");
	#pragma omp parallel for if(flag & backend::multithread)
		for( size_t n=0; n<positions.size(); ++n ) {
			const point3<T> &p (positions[n]);
			const int i (clamp<int>(p.x/dx,0,nx-2));
			const int j (clamp<int>(p.y/dx,0,ny-2));
			const int k (clamp<int>(p.z/dx,0,nz-2));
			const size_t index (k*cell_xy+j*(nx-1)+i);
			mutex[index].lock();
			hash[index].push_back(n);
			mutex[index].unlock();
		}
		printf( "Starting interpolation...\n");
		//
	#pragma omp parallel for if(flag & backend::multithread)
		for( int ij=0; ij<cell_xy; ++ij ) {
			//
			const int i = ij % (nx-1);
			const int j = ij / (nx-1);
			//
			T value_cache[6*6*6];
			std::vector<T> results_cache;
			std::vector<T> point_x_cache, point_y_cache, point_z_cache;
			//
			char head (0); int k(0);
			//
			auto incremental_load_values = [&]() {
				const int fk = clamp(k+w/2+1,0,(int)nz-1);
				for( int jj=0; jj<w; ++jj ) {
					for( int ii=0; ii<w; ++ii ) {
						const int fi = clamp(i+ii-w/2+1,0,(int)nx-1);
						const int fj = clamp(j+jj-w/2+1,0,(int)ny-1);
						value_cache[ii+jj*w+head*ww] = volume[fk*nxy+fj*nx+fi];
					}
				}
				head = (head+1) % w;
				k ++;
			};
			//
			load_initial_values(i,j,k,w,value_cache);
			//
			while( k < nz-1 ) {
				//
				const std::vector<size_t> &indices = hash[i+(nx-1)*j+cell_xy*k];
				const unsigned size (indices.size());
				results_cache.resize(size);
				point_x_cache.resize(size);
				point_y_cache.resize(size);
				point_z_cache.resize(size);
				for( unsigned n=0; n<size; ++n ) {
					point_x_cache[n] = clamp(positions[indices[n]].x/dx-i,(T)0.0,(T)1.0);
					point_y_cache[n] = clamp(positions[indices[n]].y/dx-j,(T)0.0,(T)1.0);
					point_z_cache[n] = clamp(positions[indices[n]].z/dx-k,(T)0.0,(T)1.0);
				}
				//
				switch( scheme ) {
					case filter::trilinear:
						interpolate_trilinear(
							head,size,ispc,
							value_cache,
							point_x_cache.data(),
							point_y_cache.data(),
							point_z_cache.data(),
							results_cache.data());
						break;
					case filter::WENO4:
						interpolate_WENO4(
							head,size,ispc,
							value_cache,
							point_x_cache.data(),
							point_y_cache.data(),
							point_z_cache.data(),
							results_cache.data());
						break;
					case filter::WENO6:
						interpolate_WENO6(
							head,size,ispc,
							value_cache,
							point_x_cache.data(),
							point_y_cache.data(),
							point_z_cache.data(),
							results_cache.data());
						break;
					default:
						printf("Undefined scheme! (%d)\n",(int)scheme);
						exit(1);
				}
				//
				for( unsigned n=0; n<size; ++n ) {
					results[indices[n]] = results_cache[n];
				}
				incremental_load_values();
			}
		}
	} else {
		//
	#pragma omp parallel for if(flag & backend::multithread)
		for( size_t n=0; n<positions.size(); ++n ) {
			//
			T value_cache[6*6*6];
			T result;
			T x(positions[n].x/dx), y(positions[n].y/dx), z(positions[n].z/dx);
			const int i (clamp<int>(x,0,nx-2));
			const int j (clamp<int>(y,0,ny-2));
			const int k (clamp<int>(z,0,nz-2));
			x = x-i; y = y-j; z = z-k;
			load_initial_values(i,j,k,w,value_cache);
			//
			switch( scheme ) {
				case filter::trilinear:
					interpolate_trilinear(
						0,1,ispc,value_cache,&x,&y,&z,&result);
					break;
				case filter::WENO4:
					interpolate_WENO4(
						0,1,ispc,value_cache,&x,&y,&z,&result);
					break;
				case filter::WENO6:
					interpolate_WENO6(
						0,1,ispc,value_cache,&x,&y,&z,&result);
					break;
				default:
					printf("Undefined scheme! (%d)\n",(int)scheme);
					exit(1);
			}
			results[n] = result;
		}
	}
	//
	return results;
}
//
template struct sginterp3::point3<float>;
template struct sginterp3::point3<double>;
//
template std::vector<float> sginterp3::interpolate<float>(
		const std::vector<point3<float> > &positions,
		unsigned nx, unsigned ny, unsigned nz, float dx,
		const float *volume,
		sginterp3::filter scheme,
		int flag
);
//
template std::vector<double> sginterp3::interpolate<double>(
		const std::vector<point3<double> > &positions,
		unsigned nx, unsigned ny, unsigned nz, double dx,
		const double *volume,
		sginterp3::filter scheme,
		int flag
);
//