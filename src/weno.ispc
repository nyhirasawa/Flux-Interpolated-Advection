//
static inline T sqr( T x ) {
	return x*x;
}
//
static inline uniform T sqru( uniform T x ) {
	return x*x;
}
//
static inline T weno_interp6( T x, const T v[6] ) {
	//
	uniform const T eps = 1e-30;
	//
	const T f_m2 = v[0];
	const T f_m1 = v[1];
	const T f_p0 = v[2];
	const T f_p1 = v[3];
	const T f_p2 = v[4];
	const T f_p3 = v[5];
	//
	uniform const T x_m2 = -2.0;
	uniform const T x_m1 = -1.0;
	uniform const T x_p0 = 0.0;
	uniform const T x_p1 = 1.0;
	uniform const T x_p2 = 2.0;
	uniform const T x_p3 = 3.0;
	//
	T C[3];
	C[0] = (x_p2-x) * (x_p3-x) * (1.0/20.0);
	C[1] = (x_p3-x) * (x-x_m2) * (1.0/10.0);
	C[2] = (x-x_m2) * (x-x_m1) * (1.0/20.0);
	//
	T S[3];
	S[0] = ((0814.*sqr(f_p1))+(4326.*sqr(f_p0))+(2976.*sqr(f_m1))+(0244.*sqr(f_m2))-(3579.*f_p0*f_p1)-(6927.*f_p0*f_m1)+(1854.*f_p0*f_m2)+(2634.*f_p1*f_m1)-(0683.*f_p1*f_m2)-(1659.*f_m1*f_m2)) * (1.0/180.0);
	S[1] = ((1986.*sqr(f_p1))+(1986.*sqr(f_p0))+(0244.*sqr(f_m1))+(0244.*sqr(f_p2))+(1074.*f_p0*f_p2)-(3777.*f_p0*f_p1)-(1269.*f_p0*f_m1)+(1074.*f_p1*f_m1)-(1269.*f_p2*f_p1)-(0293.*f_p2*f_m1)) * (1.0/180.0);
	S[2] = ((0814.*sqr(f_p0))+(4326.*sqr(f_p1))+(2976.*sqr(f_p2))+(0244.*sqr(f_p3))-(0683.*f_p0*f_p3)+(2634.*f_p0*f_p2)-(3579.*f_p0*f_p1)-(6927.*f_p1*f_p2)+(1854.*f_p1*f_p3)-(1659.*f_p2*f_p3)) * (1.0/180.0);
	//
	T P[3];
	P[0] = f_m2 + (f_m1-f_m2) * (x-x_m2) + (f_p0-2.*f_m1+f_m2) * (x-x_m2) * (x-x_m1) * 0.5 + (f_p1-3.*f_p0+3.*f_m1-f_m2) * (x-x_m2) * (x-x_m1) * (x-x_p0) * (1.0/6.0);
	P[1] = f_m1 + (f_p0-f_m1) * (x-x_m1) + (f_p1-2.*f_p0+f_m1) * (x-x_m1) * (x-x_p0) * 0.5 + (f_p2-3.*f_p1+3.*f_p0-f_m1) * (x-x_m1) * (x-x_p0) * (x-x_p1) * (1.0/6.0);
	P[2] = f_p0 + (f_p1-f_p0) * (x-x_p0) + (f_p2-2.*f_p1+f_p0) * (x-x_p0) * (x-x_p1) * 0.5 + (f_p3-3.*f_p2+3.*f_p1-f_p0) * (x-x_p0) * (x-x_p1) * (x-x_p2) * (1.0/6.0);
	//
	T a[3];
	T sum = 0.0;
	for( uniform int i=0; i<3; ++i ) {
		a[i] = C[i] * rcp(eps+sqr(S[i]));
		sum += a[i];
	}
	//
	return (a[0]*P[0]+a[1]*P[1]+a[2]*P[2])*rcp(sum);
};
//
static T inline weno_interp6_uniform_precompute( uniform T S[3], uniform const T v[6] ) {
	//
	uniform const T f_m2 = v[0];
	uniform const T f_m1 = v[1];
	uniform const T f_p0 = v[2];
	uniform const T f_p1 = v[3];
	uniform const T f_p2 = v[4];
	uniform const T f_p3 = v[5];
	//
	S[0] = ((0814.*sqru(f_p1))+(4326.*sqru(f_p0))+(2976.*sqru(f_m1))+(0244.*sqru(f_m2))-(3579.*f_p0*f_p1)-(6927.*f_p0*f_m1)+(1854.*f_p0*f_m2)+(2634.*f_p1*f_m1)-(0683.*f_p1*f_m2)-(1659.*f_m1*f_m2)) * (1.0/180.0);
	S[1] = ((1986.*sqru(f_p1))+(1986.*sqru(f_p0))+(0244.*sqru(f_m1))+(0244.*sqru(f_p2))+(1074.*f_p0*f_p2)-(3777.*f_p0*f_p1)-(1269.*f_p0*f_m1)+(1074.*f_p1*f_m1)-(1269.*f_p2*f_p1)-(0293.*f_p2*f_m1)) * (1.0/180.0);
	S[2] = ((0814.*sqru(f_p0))+(4326.*sqru(f_p1))+(2976.*sqru(f_p2))+(0244.*sqru(f_p3))-(0683.*f_p0*f_p3)+(2634.*f_p0*f_p2)-(3579.*f_p0*f_p1)-(6927.*f_p1*f_p2)+(1854.*f_p1*f_p3)-(1659.*f_p2*f_p3)) * (1.0/180.0);

}
static T inline weno_interp6_uniform( T x, uniform const uniform T S[3], uniform const T v[6] ) {
	//
	uniform const T eps = 1e-30;
	//
	uniform const T f_m2 = v[0];
	uniform const T f_m1 = v[1];
	uniform const T f_p0 = v[2];
	uniform const T f_p1 = v[3];
	uniform const T f_p2 = v[4];
	uniform const T f_p3 = v[5];
	//
	uniform const T x_m2 = -2.0;
	uniform const T x_m1 = -1.0;
	uniform const T x_p0 = 0.0;
	uniform const T x_p1 = 1.0;
	uniform const T x_p2 = 2.0;
	uniform const T x_p3 = 3.0;
	//
	T C[3];
	C[0] = (x_p2-x) * (x_p3-x) * (1.0/20.0);
	C[1] = (x_p3-x) * (x-x_m2) * (1.0/10.0);
	C[2] = (x-x_m2) * (x-x_m1) * (1.0/20.0);
	//
	T P[3];
	P[0] = f_m2 + (f_m1-f_m2) * (x-x_m2) + (f_p0-2.*f_m1+f_m2) * (x-x_m2) * (x-x_m1) * 0.5 + (f_p1-3.*f_p0+3.*f_m1-f_m2) * (x-x_m2) * (x-x_m1) * (x-x_p0) * (1.0/6.0);
	P[1] = f_m1 + (f_p0-f_m1) * (x-x_m1) + (f_p1-2.*f_p0+f_m1) * (x-x_m1) * (x-x_p0) * 0.5 + (f_p2-3.*f_p1+3.*f_p0-f_m1) * (x-x_m1) * (x-x_p0) * (x-x_p1) * (1.0/6.0);
	P[2] = f_p0 + (f_p1-f_p0) * (x-x_p0) + (f_p2-2.*f_p1+f_p0) * (x-x_p0) * (x-x_p1) * 0.5 + (f_p3-3.*f_p2+3.*f_p1-f_p0) * (x-x_p0) * (x-x_p1) * (x-x_p2) * (1.0/6.0);
	//
	T a[3];
	T sum = 0.0;
	for( uniform int i=0; i<3; ++i ) {
		a[i] = C[i] * rcp(eps+sqru(S[i]));
		sum += a[i];
	}
	//
	return (a[0]*P[0]+a[1]*P[1]+a[2]*P[2])*rcp(sum);
};
//
static inline T weno_interp4_uniform_precompute( uniform T S[2], uniform const T v[4] ) {
	//
	uniform const T f_m1 = v[0];
	uniform const T f_p0 = v[1];
	uniform const T f_p1 = v[2];
	uniform const T f_p2 = v[3];
	//
	S[0] = ((26.*f_p1*f_m1)-(52.*f_p0*f_m1)-(76.*f_p1*f_p0)+(25.*sqru(f_p1))+(64.*sqru(f_p0))+(13.*sqru(f_m1))) / 12.0;
	S[1] = ((26.*f_p2*f_p0)-(52.*f_p2*f_p1)-(76.*f_p1*f_p0)+(25.*sqru(f_p0))+(64.*sqru(f_p1))+(13.*sqru(f_p2))) / 12.0;
}
static inline T weno_interp4_uniform( T x, uniform const T S[2], uniform const T v[4] ) {
	//
	uniform const T eps = 1e-30;
	//
	const T f_m1 = v[0];
	const T f_p0 = v[1];
	const T f_p1 = v[2];
	const T f_p2 = v[3];
	//
	uniform const T x_m1 = -1.0;
	uniform const T x_p0 = 0.0;
	uniform const T x_p1 = 1.0;
	uniform const T x_p2 = 2.0;
	//
	T C[2];
	C[0] = (x_p2-x) * (1.0/3.0);
	C[1] = (x-x_m1) * (1.0/3.0);
	//
	T P[2];
	P[0] = f_p0 + (f_p1-f_m1) * (x-x_p0) * 0.5 + (f_p1-2.*f_p0+f_m1) * sqr(x-x_p0) * 0.5;
	P[1] = f_p0 + (-f_p2+4.*f_p1-3.*f_p0) * (x-x_p0) * 0.5 + (f_p2-2.*f_p1+f_p0) * sqr(x-x_p0) * 0.5;
	//
	T a[2];
	T sum = 0.0;
	for( uniform int i=0; i<2; ++i ) {
		a[i] = C[i] * rcp(eps+sqru(S[i]));
		sum += a[i];
	}
	//
	return (a[0]*P[0]+a[1]*P[1])*rcp(sum);
}
//
static inline T weno_interp4( T x, const T v[4] ) {
	//
	uniform const T eps = 1e-30;
	//
	const T f_m1 = v[0];
	const T f_p0 = v[1];
	const T f_p1 = v[2];
	const T f_p2 = v[3];
	//
	uniform const T x_m1 = -1.0;
	uniform const T x_p0 = 0.0;
	uniform const T x_p1 = 1.0;
	uniform const T x_p2 = 2.0;
	//
	T C[2];
	C[0] = (x_p2-x) * (1.0/3.0);
	C[1] = (x-x_m1) * (1.0/3.0);
	//
	T S[2];
	S[0] = ((26.*f_p1*f_m1)-(52.*f_p0*f_m1)-(76.*f_p1*f_p0)+(25.*sqr(f_p1))+(64.*sqr(f_p0))+(13.*sqr(f_m1))) / 12.0;
	S[1] = ((26.*f_p2*f_p0)-(52.*f_p2*f_p1)-(76.*f_p1*f_p0)+(25.*sqr(f_p0))+(64.*sqr(f_p1))+(13.*sqr(f_p2))) / 12.0;
	//
	T P[2];
	P[0] = f_p0 + (f_p1-f_m1) * (x-x_p0) * 0.5 + (f_p1-2.*f_p0+f_m1) * sqr(x-x_p0) * 0.5;
	P[1] = f_p0 + (-f_p2+4.*f_p1-3.*f_p0) * (x-x_p0) * 0.5 + (f_p2-2.*f_p1+f_p0) * sqr(x-x_p0) * 0.5;
	//
	T a[2];
	T sum = 0.0;
	for( uniform int i=0; i<2; ++i ) {
		a[i] = C[i] * rcp(eps+sqr(S[i]));
		sum += a[i];
	}
	//
	return (a[0]*P[0]+a[1]*P[1])*rcp(sum);
}
//
