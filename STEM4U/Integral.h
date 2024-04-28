// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#ifndef _STEM4U_Integral_h_
#define _STEM4U_Integral_h_

#include <ScatterDraw/DataSource.h>
#include <ScatterDraw/Equation.h>

namespace Upp {
	
enum IntegralType {TRAPEZOIDAL, SIMPSON_1_3, SIMPSON_3_8, HERMITE_3, HERMITE_5, SPLINE}; 


template <class Range>
typename Range::value_type Integral(const Range &x, const Range &y, IntegralType type = TRAPEZOIDAL) {
	using Scalar = typename Range::value_type;
	
	ASSERT(x.size() == y.size());
	if (y.size() <= 1)
		return 0;

	auto n = x.size();
	if (type == SPLINE) {
		Spline spl(x, y);
		return spl.Integral(x[0], x[n-1]);
	}
		
	Scalar ret = 0;
	
	if(type == SIMPSON_1_3 && n == 2)
		type = TRAPEZOIDAL;
	else if (type == SIMPSON_3_8) {
		if (n == 2)
			type = TRAPEZOIDAL;
		else if (n == 3)
			type = SIMPSON_1_3;
	}
	
	if (type == TRAPEZOIDAL) {
		for (int i = 1; i < n; ++i)
			ret += Avg(y[i], y[i-1])*(x[i] - x[i-1]);
	} else if (type == SIMPSON_1_3) {
		int i;
		for (i = 2; i < n; i += 2)
			ret += (x[i] - x[i-2])*(y[i-2] + 4*y[i-1] + y[i]);
		ret /= 6.;
		if (i == n)
			ret += Avg(y[n-1], y[n-2])*(x[n-1] - x[n-2]);
	} else if (type == SIMPSON_3_8) {
		int i;
		for (i = 3; i < n; i += 3)
			ret += (x[i] - x[i-3])*(y[i-3] + 3*y[i-2] + 3*y[i-1] + y[i]);
		ret /= 8.;
		if (i == n)
			ret += (x[n-1] - x[n-3])/6.*(y[n-3] + 4*y[n-2] + y[n-1]);
		else if (i == n+1)
			ret += Avg(y[n-1], y[n-2])*(x[n-1] - x[n-2]);
	} else
		NEVER();
	return ret;
}


template <class Range>
inline typename Range::value_type Calc1_3(const Range &y, typename Range::value_type dx, int n) {
	auto ret = y[0] + y[n-1];
	for (int i = 1; i < n-1; i++)
		ret += 2*y[i];
	for (int i = 1; i < n-1; i += 2)
		ret += 2*y[i];
	ret *= dx/3.;
	return ret;
}

inline double Calc1_3(Eigen::VectorXd &y, double dx, size_t n) {
	return dx/3*(y(0) + 2*(Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<2>>(y.data()+1, n/2).sum() + 
									 y.block(1, 0, n-2, 1).sum()) + y(n-1));
}

template <class Range>
typename Range::value_type Integral(const Range &y, typename Range::value_type dx, IntegralType type = TRAPEZOIDAL) {
	using Scalar = typename Range::value_type;
	
	if (y.size() <= 1)
		return 0;
	
	auto n = y.size();
	if (type == SPLINE) {
		Range x(n);
		for (int i = 0; i < n; ++i)
			x[i] = i*dx;
		Spline spl(x, y);
		return spl.Integral(x[0], x[n-1]);
	}
	
	if(type == SIMPSON_1_3 && n == 2)
		type = TRAPEZOIDAL;
	else if (type == SIMPSON_3_8) {
		if (n == 2)
			type = TRAPEZOIDAL;
		else if (n == 3)
			type = SIMPSON_1_3;
	} else if (type == HERMITE_3 && n == 2) 
		type = TRAPEZOIDAL;
	else if (type == HERMITE_5) {
		if (n == 2)
			type = TRAPEZOIDAL;
		else if (n == 3)
			type = HERMITE_3;
	}	
	
	if (type == TRAPEZOIDAL) 
		return (y.segment(1, n-2).sum() + (y(0) + y(n-1))/2)*dx;
	else if (type == SIMPSON_1_3) {
		Scalar ret0 = 0;
		if ((n-1)%2) {
			ret0 = Avg(y(n-1), y(n-2))*dx;
			--n;
		}
		return ret0 + Calc1_3(y, dx, n);
	} else if (type == SIMPSON_3_8) {
		Scalar ret0 = 0;
		int rem = (n-1)%3;
		if (rem == 2) {
			ret0 = (y[n-3] + 4*y[n-2] + y[n-1])/3.*dx;
			n -= 2;
		} else if (rem == 1) {
			ret0 = Avg(y[n-2], y[n-1])*dx;
			n--;
		}
		Scalar ret = y[0] + y[n-1];
		for (int i = 1; i < n-1; ++i)
			ret += 2*y[i];
		for (int i = 1; i < n-1; ++i) {
			ret += y[i++];
			ret += y[i++];
		}
		return ret0 + ret*dx*3./8;
	} else if (type == HERMITE_3) {
		return (y.segment(1, n-2).sum() + (y(0) + y(n-1))/2)*dx
			 + dx/24.*(3*y[0] - 4*y[1] + y[2] + y[n-3] - 4*y[n-2] + 3*y[n-1]);
	} else if (type == HERMITE_5) {
		if (n == 4)
			return ((y[0] + y[3])/2. + y[1] + y[2])*dx;
		return (y.segment(1, n-2).sum() + (y(0) + y(n-1))/2)*dx
			 - dx/144.*(25*(y[0] - y[n-1]) - 48*(y[1] - y[n-2]) + 36*(y[2] - y[n-3])
			 		  - 16*(y[3] - y[n-4]) +  3*(y[4] - y[n-5]));
	} else
		NEVER();
	return Null;
}

template <class Range>
typename Range::value_type IntegralSinCos(const Range &x, const Range &f, typename Range::value_type t, bool iscos) {
	using Scalar = typename Range::value_type;
	
	int n = f.size();
	
	ASSERT(x.size() == n);
	ASSERT(n%2);
	
	if (!(n%2))
		return Null;
	
	if (n == 0)
		return 0;
	
	Scalar dx = x[1]-x[0];
	Scalar theta = t*dx;
	Scalar sint = sin(theta);
	Scalar cost = cos(theta);
	Scalar theta2 = theta*theta;
	Scalar theta3 = theta2*theta;
	Scalar alpha, ceven, codd;
	if (6*abs(theta) <= 1) {	// Taylor series approximation
	    alpha = 2*theta3/45. - 2*PowInt(theta, 5)/315. + 2*PowInt(theta, 7)/4725.;
	    ceven =  2./3. + 2.*theta2/15. - 4.*PowInt(theta, 4)/105. + 2.*PowInt(theta, 6)/567. - 4.*PowInt(theta, 8)/22275.;
	    codd = 4./3. - 2.*theta2/15.0 + PowInt(theta, 4)/210. - PowInt(theta, 6)/11340.;
  	} else {
		alpha = (theta2 + theta*sint*cost - 2*sint*sint)/theta3;
		ceven = (2*theta + 2*theta*cost*cost - 4*sint*cost)/theta3;
		codd = 4*(sint - theta*cost)/theta3;	
	}
	Scalar even, odd;
	Scalar ret;
	if (iscos) {
		Range fcostx(n);
		for (int i = 0; i < n; ++i)
			fcostx[i] = f[i]*cos(t*x[i]);
	
		even = 0.5*f[0]*cos(t*x[0]);
		for (int i = 2; i < n-1; i += 2)
			even += fcostx[i];
		even += 0.5*fcostx[n-1];
		
		odd = 0;
		for (int i = 1; i < n-1; i += 2)
			odd += fcostx[i];
	
		ret = alpha*(f[n-1]*sin(t*x[n-1]) - f[0]*sin(t*x[0]));
	} else {
		Range fsintx(n);
		for (int i = 0; i < n; ++i)
			fsintx[i] = f[i]*sin(t*x[i]);
		
		even = 0.5*f[0]*sin(t*x[0]);
		for (int i = 2; i < n - 1; i += 2)
			even += fsintx[i];
		even += 0.5*fsintx[n-1];
		
		odd = 0;
		for (int i = 1; i < n-1; i += 2)
			odd += fsintx[i];

		ret = alpha*(f[0]*cos(t*x[0]) - f[n-1]*cos(t*x[n-1]));
	}	
	return dx*(ret + ceven*even + codd*odd);
}

}
	
#endif
