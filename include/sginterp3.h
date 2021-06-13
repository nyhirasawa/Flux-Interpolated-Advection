//
#ifndef _STREAM_GRID_INTERP3_
#define _STREAM_GRID_INTERP3_
//
#include <vector>
//
class sginterp3 {
public:
	//
	// 使用可能な補間方法
	enum filter { trilinear, WENO4, WENO6 };
	//
	// 使用可能な最適化オプション
	enum backend {
		presort = 1 << 0, // ソートしてキャッシュ効率を高める
		multithread = 1 << 1, // マルチスレッドを使用する
		simd = 1 << 2, // SIMD を利用する
	};
	//
	// 点の位置情報
	template <class T> // T は float か double
	struct point3 {
		point3() : x(0.0), y(0.0), z(0.0) {}
		point3( T x, T y, T z ) : x(x), y(y), z(z) {}
		T x; T y; T z;
	};
	//
	// 補間した値の配列が返る
	template <class T> // T は float か double
	static std::vector<T> interpolate (
		//
		// 補間点の配列
		const std::vector<point3<T> > &positions,
		//
		// ボリュームデータの解像度。ボリュームデータは nx * ny * nz の格子点を持つ。
		unsigned nx, unsigned ny, unsigned nz,
		//
		// 格子間の距離
		T dx,
		//
		// ボリュームデータ (格子点はセルの中心でなく、格子の頂点で定義されることに注意)
		// volume(i,j,k) = volume[k*(nx*ny)+j*nx+i] でアクセス可能な形式で記録されていること
		const T *volume,
		//
		// 補間方法
		filter scheme = filter::WENO6,
		//
		// 最適化オプション
		int flag = 
			backend::presort |
			backend::multithread | 
			backend::simd
	);
};
//
#endif
//