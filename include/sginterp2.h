//
#ifndef _STREAM_GRID_INTERP2_
#define _STREAM_GRID_INTERP2_
//
#include <vector>
//
class sginterp2 {
public:
	//
	// 使用可能な補間方法
	enum filter { bilinear, WENO4, WENO6 };
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
	struct point2 {
		point2() : x(0.0), y(0.0) {}
		point2( T x, T y ) : x(x), y(y) {}
		T x; T y;
	};
	//
	// 補間した値の配列が返る
	template <class T> // T は float か double
	static std::vector<T> interpolate (
		//
		// 補間点の配列
		const std::vector<point2<T> > &positions,
		//
		// ボリュームデータの解像度。ボリュームデータは nx * ny の格子点を持つ。
		unsigned nx, unsigned ny,
		//
		// 格子間の距離
		T dx,
		//
		// ボリュームデータ (格子点はセルの中心でなく、格子の頂点で定義されることに注意)
		// volume(i,j) = volume[j*nx+i] でアクセス可能な形式で記録されていること
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