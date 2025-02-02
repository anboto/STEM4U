// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#ifndef _STEM4U_utility_h_
#define _STEM4U_utility_h_

namespace Upp {

template <typename T>
bool IsNum(const Upp::Vector<T> &v) {
	if (v.IsEmpty())
		return false;
	for (const T &a : v) {
		if (!IsNum(a))
			return false;
	}
	return true;
}

template <typename T>
bool IsNum(const Upp::Array<T> &v) {
	if (v.IsEmpty())
		return false;
	for (const T &a : v) {
		if (!IsNum(a))
			return false;
	}
	return true;
}

template <class Range>
void CleanOutliers(const Range &x, const Range &y, const Range &filtery, Range &rretx, Range &rrety, 
				   const typename Range::value_type& ratio, const typename Range::value_type& zero = 0) {
	ASSERT(x.size() == y.size() && x.size() == filtery.size());
 
 	Range retx(x.size()), rety(x.size());
 
 	int id = 0;
	for (int i = 0; i < x.size(); ++i) {
		if (EqualRatio(y[i], filtery[i], ratio, zero)) {
			retx[id] = x[i];
			rety[id] = y[i];
			id++;
		}
	}
	ResizeConservative(retx, id); 
	ResizeConservative(rety, id);
	
	rretx = pick(retx);
	rrety = pick(rety);
}

template <class Range>
void CleanCondition(const Range &x, const Range &y, Range &rretx, Range &rrety, 
		Function <bool(int id)> Cond) {
	ASSERT(x.size() == y.size());
 
 	Range retx(x.size()), rety(x.size());
 	
 	int id = 0;
	for (int i = 0; i < x.size(); ++i) {
		if (Cond(i)) {
			retx[id] = x[i];
			rety[id] = y[i];
			id++;
		}
	}
	ResizeConservative(retx, id); 
	ResizeConservative(rety, id);

	rretx = pick(retx);
	rrety = pick(rety);
}

template <class Range>	// Gets the most probable sample rate, or the average if the most probable probability is lower than rate
typename Range::value_type GetSampleRate(const Range &x, int numDecimals, double rate) {
	using Scalar = typename Range::value_type;
	
	int n = int(x.size());
	if (n < 2)
		return Null;
	
	Vector<Scalar> delta;
	Vector<int> num;
	
	for (int i = 1; i < n; ++i) {
		Scalar d = x[i]-x[i-1];
		int id = FindRoundDecimals(delta, d, numDecimals);	
		if (id >= 0)
			num[id]++;
		else {
			delta << d;
			num << 1;		
		}
	}
	int nummx = num[0], idmx = 0;
	for (int i = 1; i < delta.size(); ++i) {
		if (num[i] > nummx) {
			nummx = num[i];
			idmx = i;
		}
	}
	if (num[idmx]/double(n-1) > rate)
		return delta[idmx];
	else {
		Scalar avg = 0;
		for (int i = 1; i < n; ++i) 
			avg += (x[i]-x[i-1]);
		return avg/(n-1);
	}
}

template <typename T>
void LinSpaced(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, int n, T min, T max) {
	v.setLinSpaced(n, min, max);
}

template <class Range>
void LinSpaced(Range &v, int n, typename Range::value_type min, typename Range::value_type max) {
	ASSERT(n > 0);
	Resize(v, n);
	if (n == 1)
		v[0] = min;
	else {
		typename Range::value_type d = (max - min)/(n - 1);
		for (int i = 0; i < n; ++i)
			v[i] = min + d*i;
	}
}

template <typename T>
void Arange(Eigen::Matrix<T, Eigen::Dynamic, 1> &v, T min, T max, T step = 1) {
	int num = int((max - min)/step)+1;
	v = Eigen::Matrix<T, Eigen::Dynamic, 1>::LinSpaced(num, min, min + (num-1)*step);
}

template <class Range>
void Arange(Range &v, typename Range::value_type min, typename Range::value_type max, typename Range::value_type step = 1) {
	int num = int((max - min)/step)+1;
	LinSpaced(v, num, min, min + (num-1)*step);
}


template <class Range>
void CircShift(const Range& in, int len, Range &out) {
	std::rotate_copy(in, in + len, End(in), out);
}

template <typename T>
void CircShift(const Eigen::Matrix<T, Eigen::Dynamic, 1> &in, int len, Eigen::Matrix<T, Eigen::Dynamic, 1> &out) {
	Resize(out, in.size());
	out.segment(len, in.size() - len) = in.segment(0, in.size() - len); 
	out.segment(0, len) = in.segment(in.size() - len, len);
}

template <class Range>
void NextPow2(const Range& in, Range &out) {
	Resize(out, in.size());
	for (int i = 0; i < in.size(); ++i)
		out[i] = ceil(log(abs(in[i])))/log(2);
}

template <typename T>
T NextPow2(const T& in) {
	ASSERT(in > 0);
	return T(ceil(log(in))/log(2));
}

template <typename T>
inline T Avg(const Eigen::Matrix<T, Eigen::Dynamic, 1> &d) {
	if (d.size() == 0)
		return Null;
	return d.mean();
}

template <class Range>
inline typename Range::value_type Avg(const Range &d) {
	if (d.size() == 0)
		return Null;
	return std::accumulate(d.begin(), d.end(), 0.)/d.size();
}

template <class Range>
typename Range::value_type R2(const Range &serie, const Range &serie0, typename Range::value_type meanserie = Null) {
	using Scalar = typename Range::value_type;
	
	if (IsNull(meanserie))
		meanserie = Avg(serie);
	Scalar sse = 0, sst = 0;
	auto sz = min(serie.size(), serie0.size());
	for (auto i = 0; i < sz; ++i) {
		auto y = serie[i];
		auto err = y - serie0[i];
		sse += err*err;
		auto d = y - meanserie;
		sst += d*d;
	}
	if (sst < 1E-50 || sse > sst)
		return 0;
	return 1 - sse/sst;
}

template <class Range>
typename Range::value_type R2(const Range &tserie, const Range &serie, const Range &tserie0, const Range &serie0, typename Range::value_type meanserie = Null) {
	Range nserie0;
	Resample(tserie0, serie0, tserie, nserie0);
	
	return R2(serie, nserie0, meanserie);
}

template <typename T>
inline T RMS(const Eigen::Matrix<T, Eigen::Dynamic, 1> &d) {
	return sqrt(d.array().square().sum()/d.size());
}

template <class Range>
typename Range::value_type RMS(const Range &d) {
	using Scalar = typename Range::value_type;
	
	Scalar accum = 0;
	std::for_each (std::begin(d), std::end(d), [&](const Scalar v) {
    	accum += sqr(v);
	});
	return sqrt(accum/d.size());
}

template <class Range>
typename Range::value_type RMSE(const Range &serie, const Range &serie0) {
	ASSERT(serie.size() == 0 || serie.size() == serie0.size());
	using Scalar = typename Range::value_type;
	
	Scalar ret = 0;
	auto sz = min(serie.size(), serie0.size());
	for (auto i = 0; i < sz; ++i) 
		ret += sqr(serie[i] - serie0[i]);

	if (sz == 0)
		return Null;
	return sqrt(ret/sz);
}

template <class Range>
typename Range::value_type RMSE(const Range &tserie, const Range &serie, const Range &tserie0, const Range &serie0) {
	Range nserie0;
	ResampleY(tserie0, serie0, tserie, nserie0);
	
	return RMSE(serie, nserie0);
}

template <class Range>
typename Range::value_type MAE(const Range &serie, const Range &serie0) {
	ASSERT(serie.size() == 0 || serie.size() == serie0.size());
	using Scalar = typename Range::value_type;
	
	Scalar ret = 0;
	auto sz = min(serie.size(), serie0.size());
	for (auto i = 0; i < sz; ++i) 
		ret += abs(serie[i] - serie0[i]);

	if (sz == 0)
		return Null;
	return ret/sz;
}

template <class Range>
typename Range::value_type MAE(const Range &tserie, const Range &serie, const Range &tserie0, const Range &serie0) {
	Range nserie0;
	ResampleY(tserie0, serie0, tserie, nserie0);
	
	return MAE(serie, nserie0);
}

template <typename T>
inline T StdDev(const Eigen::Matrix<T, Eigen::Dynamic, 1> &d) {
	if (d.size() < 2)
		return Null;
	return sqrt((d.array() - d.mean()).square().sum()/(d.size() - 1));
}

template <class Range>
inline typename Range::value_type StdDev(const Range &d) {
	if (d.size() < 2)
		return Null;
	
	using Scalar = typename Range::value_type;
	
	Scalar avg = Avg(d);

	Scalar accum = 0;
	std::for_each (std::begin(d), std::end(d), [&](const Scalar v) {
    	accum += sqr(v - avg);
	});
	return sqrt(accum/(d.size() - 1));
}

template <class Range>
void Segment(const Range &x, const Range &y, typename Range::value_type fromx, typename Range::value_type tox, Range &xx, Range &yy) {
	ASSERT(x.size() == y.size());
	ASSERT(fromx <= tox);
	
	int ibegin = Null;
	for (int i = 0; i < x.size(); ++i) {
		if (x[i] >= fromx) {
			ibegin = i;
			break;
		}
	}
	if (IsNull(ibegin)) {
		Clear(xx);
		Clear(yy);
	}
	
	int iend = Null;
	for (int i = int(x.size())-1; i >= 0; --i) {
		if (x[i] <= tox) {
			iend = i;
			break;
		}
	}
	if (IsNull(iend)) {
		Clear(xx);
		Clear(yy);
	}
	
	xx = Segment(x, ibegin, iend-ibegin+1);
	yy = Segment(y, ibegin, iend-ibegin+1);
}

template <class Range>
typename Range::value_type SawTeethRatio(const Range &d) {
	using Scalar = typename Range::value_type;
	
	if (d.size() < 3)
		return 0;
	
	Scalar mean = Avg(d);
	int numcrosses = 0;
	for (int i = 1; i < d.size(); ++i) {
		if ((d[i] > mean && d[i-1] <= mean) || (d[i] < mean && d[i-1] >= mean)) 
			numcrosses++;
	}
	return Scalar(numcrosses)/(d.size()-1);
}


// y = ax + b
template <class Range>
void LinearFitting(const Range &x, const Range &y, typename Range::value_type &a, typename Range::value_type &b, const int clen = Null) {
	int len = IsNull(clen) ? x.size() : clen;
	ASSERT(len <= x.size() && len <= y.size());
	using Scalar = typename Range::value_type;

	a = 0;	
	Scalar sumx = 0, sumy = 0;	
	for (int i = 0; i < len; i++) { 
		sumx += x[i];
		sumy += y[i];
	}
	Scalar sum = len;
	Scalar sxoss = sumx/sum;
	Scalar sumtotal = 0;
	for (int i = 0; i < len; i++) {
		Scalar ti = x[i] - sxoss;
		sumtotal += ti*ti;
		a += ti*y[i];
	}
	a /= sumtotal; 
	b = (sumy - sumx*a)/sum;
}

#define LinearRegression	LinearFitting

// https://en.wikipedia.org/wiki/Smoothstep
template <typename T>
T SmoothStep(T x, int order = 3, T x0 = 0, T x1 = 1, T y0 = 0, T y1 = 1) {
	ASSERT (order == 3 || order == 5 || order == 7);
	if (x < x0)
		return y0;
	if (x > x1)
		return y1;
	
	double r = double(x - x0)/double(x1 - x0);
	
	T ret = Null;
	switch (order) {
	case 3: 	ret = T(r*r*(3 - 2*r));								break;
	case 5: 	ret = T(r*r*r*(10 - 15*r + 6*r*r));					break;
	case 7:		ret = T(r*r*r*r*(35 - 84*r + 70*r*r - 20*r*r*r));	break;
	default:	NEVER();
	}
	return ret*(y1 - y0) + y0;
}

template <typename In>
class Homography {
public:
	Homography() {}

	void QuadToQuad(const Point_<In> &from0, const Point_<In> &from1, const Point_<In> &from2, const Point_<In> &from3,
			   		const Point_<In> &to0,   const Point_<In> &to1,   const Point_<In> &to2,   const Point_<In> &to3) {
	    Eigen::Matrix<double, 4, 1> x, y, u, v;
	    
	    x << from0.x, from1.x, from2.x, from3.x;
	    y << from0.y, from1.y, from2.y, from3.y;
	    u << to0.x,   to1.x,   to2.x,   to3.x;
	    v << to0.y,   to1.y,   to2.y,   to3.y;

    	Eigen::Matrix<double, 9, 9> A = Eigen::Matrix<double, 9, 9>::Zero();
	    int j = 0;
	
	    for (int i = 0; i < 4; ++i) {
	        A(j, 0) = -x(i); 			A(j, 1) = -y(i); 			A(j, 2) = -1;
	        A(j, 6) = x(i) * u(i); 		A(j, 7) = y(i) * u(i); 		A(j, 8) = u(i);
	
	        A(j + 1, 3) = -x(i); 		A(j + 1, 4) = -y(i); 		A(j + 1, 5) = -1;
	        A(j + 1, 6) = x(i) * v(i); 	A(j + 1, 7) = y(i) * v(i); 	A(j + 1, 8) = v(i);
	
	        j += 2;
	    }
	    A(8, 8) = 1;  // Assuming h_9 = 1

	    Eigen::Matrix<double, 9, 1> b;
	    b << 0, 0, 0, 0, 0, 0, 0, 0, 1;

    	Eigen::Matrix<double, 9, 1> h = A.colPivHouseholderQr().solve(b);
    
    	H = Eigen::Map<Eigen::Matrix<double, 3, 3>>(h.data()).transpose();
	}

	void QuadToRectangle(const Point_<In> &from0, const Point_<In> &from1, const Point_<In> &from2, const Point_<In> &from3,
			   			 const Size_<In> &sz) {
		QuadToQuad(from0, from1, from2, from3, Point_<In>(0, 0), Point_<In>(sz.cx, 0), Point_<In>(sz.cx, sz.cy), Point_<In>(0, sz.cy));
	}
			   			
	Pointf Transform(const Point_<In> &from) {
		Eigen::Vector3d p;
		p << from.x, from.y, 1;
		
		Eigen::Vector2d p1 = Eigen::Vector3d(H * p).hnormalized();
		return Point_<double>(p1(0), p1(1));
	}
	
private:
	Eigen::Matrix3d H;
};

Image ApplyHomography(const Image& orig, const Color &back, 
			const Point &from0, const Point &from1, const Point &from2, const Point &from3,
			const Point &to0,   const Point &to1,   const Point &to2,   const Point &to3, bool bilinear = true, bool fit = true);

Image ApplyHomography(const Image& orig, const Color &back, 
			const Point &from0, const Point &from1, const Point &from2, const Point &from3, const Size &sz, bool bilinear = true);

RGBA GetPixelBilinear(const Image &img, double x, double y);
RGBA GetPixel(const Image &img, double x, double y, bool bilinear);
			
}

#endif
