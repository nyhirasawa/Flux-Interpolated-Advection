/*
**	facecutter.h
**
**	Created by Ryoichi Ando <rand@nii.ac.jp> on December 31, 2020.
**
**	Permission is hereby granted, free of charge, to any person obtaining a copy of
**	this software and associated documentation files (the "Software"), to deal in
**	the Software without restriction, including without limitation the rights to use,
**	copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
**	Software, and to permit persons to whom the Software is furnished to do so,
**	subject to the following conditions:
**
**	The above copyright notice and this permission notice shall be included in all copies
**	or substantial portions of the Software.
**
**	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
**	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
**	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
**	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
**	CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
**	OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <functional>
#include <limits>
#include <vector>
#include <random>
//
#define STACKFREE		0
#define MAX_TMP_POINTS	32
#define MAX_STACK_PNUM	32
#define MAX_STACK_VNUM	128
//
struct point2 {
	float x[2];
	point2() { x[0] = x[1] = 0.0; }
	point2 ( float x, float y ) {
		this->x[0] = x; this->x[1] = y;
	}
	point2 operator*(float s) const { return point2(s*x[0],s*x[1]); }
	point2 operator+(const point2 &a) const { return point2(x[0]+a.x[0],x[1]+a.x[1]); }
	point2 operator-(const point2 &a) const { return point2(x[0]-a.x[0],x[1]-a.x[1]); }
};
//
static inline point2 operator*(float s, const point2 &p) {
	return p*s;
}
//
struct vertex2 {
	vertex2 () {}
	vertex2 ( const point2 &p, const point2 &c ) { this->p = p; this->c = c; }
	point2 p;	// 頂点の座標
	point2 c;	// 頂点カラー
};
//
static inline float sgn( float val ) {
	return val < 0.0 ? -1.0 : 1.0;
}
//
static void cut_face( 
	const vertex2 *_polygon,		// 入力のポリゴンの頂点の先頭のポインタ
	const unsigned char _num,		// 入力のポリゴンの頂点の数
	std::function<void(const vertex2 *polygon, unsigned char num)> out_func, // ポリゴン出力関数
	const unsigned char _dim=0,
	const bool debug=false
) {
	//
#if STACKFREE
	vertex2 vertices_stack[MAX_STACK_VNUM];
	unsigned char polygon_num_stack[MAX_STACK_PNUM];
	unsigned char dim_stack[MAX_STACK_PNUM];
	unsigned vertex_head (0);
	unsigned polygon_head (0);
	//
	auto push_polygon = [&]( const vertex2 *polygon, const unsigned char num, const unsigned char dim ) {
		//
		assert( vertex_head+num < MAX_STACK_VNUM );
		assert( polygon_head+1 < MAX_STACK_PNUM );
		//
		std::memcpy( vertices_stack+vertex_head, polygon, num*sizeof(vertex2) );
		polygon_num_stack[polygon_head] = num;
		dim_stack[polygon_head] = dim;
		//
		vertex_head += num;
		polygon_head ++;
	};
	//
	auto pop_polygon = [&]( vertex2 *out_polygon, unsigned char *out_num, unsigned char *out_dim ) {
		//
		if( polygon_head ) {
			polygon_head --;
			*out_num = polygon_num_stack[polygon_head];
			*out_dim = dim_stack[polygon_head];
			vertex_head -= *out_num;
			std::memcpy(out_polygon,vertices_stack+vertex_head,(*out_num)*sizeof(vertex2) );
			return true;
		} else {
			return false;
		}
	};
	//
	push_polygon(_polygon,_num,_dim);
#endif
	//
#if STACKFREE
	vertex2 polygon[MAX_TMP_POINTS];
	unsigned char num, dim;
	while( pop_polygon(polygon,&num,&dim) ) {
#else
	const vertex2 *polygon (_polygon);
	const unsigned char num (_num);
	const unsigned char dim (_dim);
	{
#endif
		//
		static const float float_max = std::numeric_limits<float>::max();
		static const float float_lowest = std::numeric_limits<float>::lowest();
		//
		float min_x = std::numeric_limits<float>::max();
		float max_x = std::numeric_limits<float>::lowest();
		for( int n=0; n<num; ++n ) {
			min_x = std::min(min_x,polygon[n].c.x[dim]);
			max_x = std::max(max_x,polygon[n].c.x[dim]);
			if( debug ) {
				printf( "v[%d] = (%g,%g) (%g)\n", n, polygon[n].p.x[0], polygon[n].p.x[1], polygon[n].c.x[dim] );
			}
		}
		//
		int start = std::ceil(min_x);
		int end = std::floor(max_x);
		if( min_x - start == 0.0 ) start += 1;
		if( max_x - end == 0.0 ) end -= 1;
		//
		bool cut_produced (false);
		if( debug ) printf( "min_x = %.2f, max_x = %.2f, start=%d, end=%d\n", min_x, max_x, start, end );
		if( start <= end ) {
			//
			const float iso_value (start);
			if( debug ) printf( "iso_value = %.2f\n", iso_value );
			//
			vertex2 vertex_array[2][MAX_TMP_POINTS];
			int vertex_array_head[2] = { 0, 0 };
			int vertex_array_slot (0);
			//
			vertex_array[0][vertex_array_head[0]++] = polygon[0];
			for( int n=0; n<num; ++n ) {
				//
				const int m = (n+1) % num;
				const float v0 = polygon[n].c.x[dim] - iso_value;
				const float v1 = polygon[m].c.x[dim] - iso_value;
				//
				if( (v0 || v1) and sgn(v0) * sgn(v1) < 0.0 ) {
					const float inv_len(1.0/(v1-v0));
					const float t0 (v1 * inv_len), t1 (-v0 * inv_len);
					cut_produced = true;
					//
					vertex2 vertex_add(
						t0 * polygon[n].p + t1 * polygon[m].p,
						point2()
					);
					vertex_add.c.x[dim] = iso_value;
					vertex_add.c.x[1-dim] = t0 * polygon[n].c.x[1-dim] + t1 * polygon[m].c.x[1-dim];
					//
					if( vertex_array_slot == 0 ) {
						vertex_array[0][vertex_array_head[0]++] = vertex_add;
						vertex_array_slot = 1;
						vertex_array_head[1] = 0;
						if( t1 < 1.0 ) {
							vertex_array[vertex_array_slot][vertex_array_head[vertex_array_slot]++] = vertex_add;
						}
					} else {
						vertex_array[1][vertex_array_head[1]++] = vertex_add;
#if STACKFREE
						push_polygon(vertex_array[1],vertex_array_head[1],dim);
#else
						cut_face(vertex_array[1],vertex_array_head[1],out_func,dim);
#endif
						vertex_array_slot = 0;
						if( t1 < 1.0 ) {
							vertex_array[vertex_array_slot][vertex_array_head[vertex_array_slot]++] = vertex_add;
						}
					}
				}
				//
				if( n < num-1 ) {
					vertex_array[vertex_array_slot][vertex_array_head[vertex_array_slot]++] = polygon[m];
				}
				//
				assert( vertex_array_head[0] < MAX_TMP_POINTS );
				assert( vertex_array_head[1] < MAX_TMP_POINTS );
			}
			if( cut_produced ) {
				if( vertex_array_slot != 0 ) {
					int cut_count (0);
					for( int n=0; n<num; ++n ) {
						const int m = (n+1) % num;
						const float v0 = polygon[n].c.x[dim] - iso_value;
						const float v1 = polygon[m].c.x[dim] - iso_value;
						if( (v0 || v1) and sgn(v0) * sgn(v1) < 0.0 ) {
							++ cut_count;
						}
						printf( "n = %d, v0 = %f, v1 = %f, cut_count = %d\n", n, v0, v1, cut_count );
					}
					exit(0);
				}
#if STACKFREE
				push_polygon(vertex_array[0],vertex_array_head[0],dim);
#else
				cut_face(vertex_array[0],vertex_array_head[0],out_func,dim);
#endif
			}
		}
		if( ! cut_produced ) {
			if( dim == 0 ) {
#if STACKFREE
				push_polygon(polygon,num,1);
#else
				cut_face(polygon,num,out_func,1);
#endif
			} else {
				out_func(polygon,num);
			}
		}
	}
}
//
static void get_center_of_gravity(
	const point2 *points,
	const unsigned char num,
	point2 &center, float &area ) {
	//
	center = point2();
	area = 0.0;
	for( unsigned char n=0; n<num; ++n ) {
		//
		const point2 &p0 = points[n];
		const point2 &p1 = points[(n+1)%num];
		const float a = p0.x[0]*p1.x[1]-p1.x[0]*p0.x[1];
		area += a;
		center.x[0] += a*(p0.x[0]+p1.x[0]);
		center.x[1] += a*(p0.x[1]+p1.x[1]);
	}
	//
	area = std::abs(0.5 * area);
	if( area ) {
		center = center * (1.0 / (6.0 * area));
	} else {
		center = point2();
		for( unsigned char n=0; n<num; ++n ) {
			center = center + points[n];
		}
		center = center * (1.0 / num);
	}
}
//
static void subdivide(
	const point2 *polygon,
	const unsigned char num,
	const float minimal_area,
	std::function<void(
		const point2 *polygon,
		unsigned char num,
		const point2 &center,
		const float area )> output_func
) {
	//
	auto sqr = []( float x ) { return x*x; };
	//
	// 再帰的に三角形を分割する関数
	std::function<void(const point2 *points)> recursive_triangle_split = 
	[&](const point2 *points) {
		//
		point2 center; float area;
		get_center_of_gravity(points,3,center,area);
		//
		if( area > minimal_area ) {
			//
			unsigned char max_index (0);
			float min_len = std::numeric_limits<float>::max();
			float max_len = std::numeric_limits<float>::lowest();
			//
			for( unsigned char n=0; n<3; ++n ) {
				const point2 &p0 = points[n];
				const point2 &p1 = points[(n+1)%3];
				const float len = sqr(p0.x[0]-p1.x[0])+sqr(p0.x[1]-p1.x[1]);
				min_len = std::min(min_len,len);
				if( max_len < len ) {
					max_len = len;
					max_index = n;
				}
			}
			//
			const point2 &p0 = points[max_index];
			const point2 &p1 = points[(max_index+1)%3];
			const point2 &p2 = points[(max_index+2)%3];
			//
			const point2 tri0[] = { p0, 0.5*(p0+p1), p2 };
			const point2 tri1[] = { 0.5*(p0+p1), p1, p2 };
			//
			recursive_triangle_split(tri0);
			recursive_triangle_split(tri1);
			//
		} else if( area ) {
			output_func(points,3,center,area);
		}
	};
	//
	// 多角形を三角形ポリゴンに分割する関数
	if( num == 3 ) {
		recursive_triangle_split(polygon);
	} else {
		//
		point2 center; float area;
		get_center_of_gravity(polygon,num,center,area);
		//
		if( area > minimal_area ) {
			for( unsigned char n=0; n<num; ++n ) {
				const point2 tri[] = { center, polygon[n], polygon[(n+1)%num] };
				recursive_triangle_split(tri);
			}
		} else if( area ) {
			output_func(polygon,num,center,area);
		}
	}
}