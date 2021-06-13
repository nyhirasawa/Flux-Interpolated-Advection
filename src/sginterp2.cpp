//
#include "sginterp2.h"
#include "weno.h"
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <algorithm>
//
extern "C" {
	//
	void ispc_interpolate_bilinear_float2(
		unsigned char head, unsigned size,
		const float value_cache[],
		const float x_array[],
		const float y_array[],
		float result[] );
	//
	void ispc_interpolate_bilinear_double2(
		unsigned char head, unsigned size,
		const double value_cache[],
		const double x_array[],
		const double y_array[],
		double result[] );
	//
	void ispc_interpolate_weno4_float2(
		unsigned char head, unsigned size,
		const float value_cache[],
		const float x_array[],
		const float y_array[],
		float result[] );
	//
	void ispc_interpolate_weno4_double2(
		unsigned char head, unsigned size,
		const double value_cache[],
		const double x_array[],
		const double y_array[],
		double result[] );
	//
	void ispc_interpolate_weno6_float2(
		unsigned char head, unsigned size,
		const float value_cache[],
		const float x_array[],
		const float y_array[],
		float result[] );
	//
	void ispc_interpolate_weno6_double2(
		unsigned char head, unsigned size,
		const double value_cache[],
		const double x_array[],
		const double y_array[],
		double result[] );
}
//
template <class T>
static T clamp( T v, T min, T max ) {
	return std::max(min,std::min(v,max));
}
//
template <class T>
static void interpolate_bilinear(
	unsigned char head0, unsigned size, bool ispc,
	const T *value_cache,
	const T *x_array,
	const T *y_array,
	T *result ) {
	//
	if( ispc ) {
		if( std::is_same<T,float>::value ) {
			ispc_interpolate_bilinear_float2(head0,size,
				(const float *)value_cache,
				(const float *)x_array,(const float *)y_array,
				(float *)result);
		} else if( std::is_same<T,double>::value ) {
			ispc_interpolate_bilinear_double2(head0,size,
				(const double *)value_cache,
				(const double *)x_array,(const double *)y_array,
				(double *)result);
		}
	} else {
		//
		const unsigned char w (2);
		const unsigned char head1 ((head0+1)%w);
		for( unsigned n=0; n<size; ++n ) {
			//
			const T x (x_array[n]), y (y_array[n]);
			//
			const T &v00 = value_cache[w*head0+0];
			const T &v10 = value_cache[w*head0+1];
			const T &v01 = value_cache[w*head1+0];
			const T &v11 = value_cache[w*head1+1];
			//
			result[n] = 
				(1.0-x)*(1.0-y)*v00 + (x)*(1.0-y)*v10 +
				(1.0-x)*(y)*v01 + (x)*(y)*v11;
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
	T *result ) {
	//
	if( ispc ) {
		if( std::is_same<T,float>::value ) {
			ispc_interpolate_weno4_float2(head0,size,
				(const float *)value_cache,
				(const float *)x_array,(const float *)y_array,
				(float *)result);
		} else if( std::is_same<T,double>::value ) {
			ispc_interpolate_weno4_double2(head0,size,
				(const double *)value_cache,
				(const double *)x_array,(const double *)y_array,
				(double *)result);
		}
	} else {
		//
		const unsigned char w (4);
		for( unsigned n=0; n<size; ++n ) {
			//
			const T x (x_array[n]), y (y_array[n]);
			//
			T tmp_x[4];
			for( unsigned m=0; m<4; ++m ) {
				const unsigned char head ((head0+m)%w);
				tmp_x[m] = WENO<T>::interp4(x,value_cache+head*w);
			}
			result[n] = WENO<T>::interp4(y,tmp_x);
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
	T *result ) {
	//
	if( ispc ) {
		if( std::is_same<T,float>::value ) {
			ispc_interpolate_weno6_float2(head0,size,
				(const float *)value_cache,
				(const float *)x_array,(const float *)y_array,
				(float *)result);
		} else if( std::is_same<T,double>::value ) {
			ispc_interpolate_weno6_double2(head0,size,
				(const double *)value_cache,
				(const double *)x_array,(const double *)y_array,
				(double *)result);
		}
	} else {
		//
		const unsigned char w (6);
		for( unsigned n=0; n<size; ++n ) {
			//
			const T x (x_array[n]), y (y_array[n]);
			//
			T tmp_x[6];
			for( unsigned m=0; m<6; ++m ) {
				const unsigned char head ((head0+m)%w);
				tmp_x[m] = WENO<T>::interp6(x,value_cache+head*w);
			}
			result[n] = WENO<T>::interp6(y,tmp_x);
		}
	}
}
//
template <class T>
std::vector<T> sginterp2::interpolate(
		const std::vector<point2<T> > &positions,
		unsigned nx, unsigned ny, T dx,
		const T *volume,
		filter scheme,
		int flag
) {
	//
	std::vector<T> results(positions.size());
	//
	const bool ispc = flag & backend::simd;
	unsigned char w (0);
	switch( scheme ) {
		case filter::bilinear:
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
	//
	auto load_initial_values = [&]( int i, int j, int w, T *value_cache ) {
		for( int ii=0; ii<w; ++ii ) {
			for( int jj=0; jj<w; ++jj ) {
				const int fi = clamp(i+ii-w/2+1,0,(int)nx-1);
				const int fj = clamp(j+jj-w/2+1,0,(int)ny-1);
				value_cache[ii+jj*w] = volume[fj*nx+fi];
			}
		}
	};
	//
	if( flag & backend::presort ) {
		//
		std::vector<std::vector<size_t> > hash((nx-1)*(ny-1));
		std::vector<std::mutex> mutex((nx-1)*(ny-1));
		//
	#pragma omp parallel for if(flag & backend::multithread)
		for( size_t n=0; n<positions.size(); ++n ) {
			const point2<T> &p (positions[n]);
			const int i (clamp<int>(p.x/dx,0,nx-2));
			const int j (clamp<int>(p.y/dx,0,ny-2));
			const size_t index(j*(nx-1)+i);
			mutex[index].lock();
			hash[index].push_back(n);
			mutex[index].unlock();
		}
		//
	#pragma omp parallel for if(flag & backend::multithread)
		for( int i=0; i<nx-1; ++i ) {
			//
			T value_cache[6*6];
			std::vector<T> results_cache, point_x_cache, point_y_cache;
			char head (0); int j(0);
			//
			auto incremental_load_values = [&]() {
				const int fj = clamp(j+w/2+1,0,(int)ny-1);
				for( int ii=0; ii<w; ++ii ) {
					const int fi = clamp(i+ii-w/2+1,0,(int)nx-1);
					value_cache[ii+head*w] = volume[fj*nx+fi];
				}
				head = (head+1) % w;
				j ++;
			};
			//
			load_initial_values(i,j,w,value_cache);
			//
			while( j < ny-1 ) {
				//
				const std::vector<size_t> &indices = hash[i+(nx-1)*j];
				const unsigned size (indices.size());
				results_cache.resize(size);
				point_x_cache.resize(size);
				point_y_cache.resize(size);
				for( unsigned n=0; n<size; ++n ) {
					point_x_cache[n] = clamp(positions[indices[n]].x/dx-i,(T)0.0,(T)1.0);
					point_y_cache[n] = clamp(positions[indices[n]].y/dx-j,(T)0.0,(T)1.0);
				}
				//
				switch( scheme ) {
					case filter::bilinear:
						interpolate_bilinear(
							head,size,ispc,
							value_cache,
							point_x_cache.data(),
							point_y_cache.data(),
							results_cache.data());
						break;
					case filter::WENO4:
						interpolate_WENO4(
							head,size,ispc,
							value_cache,
							point_x_cache.data(),
							point_y_cache.data(),
							results_cache.data());
						break;
					case filter::WENO6:
						interpolate_WENO6(
							head,size,ispc,
							value_cache,
							point_x_cache.data(),
							point_y_cache.data(),
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
			T value_cache[6*6];
			T result, x(positions[n].x/dx), y(positions[n].y/dx);
			const int i (clamp<int>(x,0,nx-2));
			const int j (clamp<int>(y,0,ny-2));
			//
			x = x-i; y = y-j;
			load_initial_values(i,j,w,value_cache);
			//
			switch( scheme ) {
				case filter::bilinear:
					interpolate_bilinear(
						0,1,ispc,value_cache,&x,&y,&result);
					break;
				case filter::WENO4:
					interpolate_WENO4(
						0,1,ispc,value_cache,&x,&y,&result);
					break;
				case filter::WENO6:
					interpolate_WENO6(
						0,1,ispc,value_cache,&x,&y,&result);
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
template struct sginterp2::point2<float>;
template struct sginterp2::point2<double>;
//
template std::vector<float> sginterp2::interpolate<float>(
		const std::vector<point2<float> > &positions,
		unsigned nx, unsigned ny, float dx,
		const float *volume,
		sginterp2::filter scheme,
		int flag
);
//
template std::vector<double> sginterp2::interpolate<double>(
		const std::vector<point2<double> > &positions,
		unsigned nx, unsigned ny, double dx,
		const double *volume,
		sginterp2::filter scheme,
		int flag
);
//