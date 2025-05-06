// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2025, the Anboto author and contributors
#include <Core/Core.h>
#include <Functions4U/Functions4U.h>

using namespace Upp;


#include <STEM4U/IntInf.h>
#include <STEM4U/Rational.h>
#include <STEM4U/Polynomial.h>


void TestPolynomial() {
	UppLog() << "\n\nRational demo. A rational number based on an arbitrary precision integer";
	
	int n = 6;
	Rational NT = 111;
	
	int r = int(Pow10Int<double>(int(log10(int(NT))-1)));
	int m = int(NT/2.);

	int M = int((NT-1.) / 2.);

	double csi_n_num = 2*(2*n+1);
	double csi_n_den = n+1;
	
	Rational gamma_n_num = NT;
	double gamma_n_den = 2 * n + 1;
	
	for (int j = 1; j < n+1; ++j) 
		gamma_n_num *= NT*NT - double(j*j);
	
	Rational gamma_n = gamma_n_num/gamma_n_den;
	
	Vector<Polynomial<Rational>> q;
	
	q << Polynomial<Rational>(1);
	UppLog() << "\n" << q[0];
	
	q << Polynomial<Rational>(0, 2);
	UppLog() << "\n" << q[1];

	for (int i = 2; i < n+2; ++i) {
		Rational ii = i;
		
		q << Polynomial<Rational>(0, (2.*ii - 1.)*2./ii) * q[i-1] - q[i-2] * (((ii-1.)*(NT*NT - (ii-1.)*(ii-1.)))/double(i));
		UppLog() << "\n" << q[i];
	}
	
	auto sg = q[n+1].Order(-1);
	UppLog() << "\n" << "sg " << sg;
	
	auto dsg = sg.Diff();	
	UppLog() << "\n" << "dsg " << dsg;

	auto num = q[n].y(0).Simplify(); 
	auto den = ((gamma_n * csi_n_num) / csi_n_den).Simplify();
	UppLog() << "\n" << "num: " << num;
	UppLog() << "\n" << "den: " << den;

	Vector<Rational> b;	
	b.SetCount(int(NT), 0);
 
	Rational sum_bN = 0;
	for (int l = -m; l < 1; ++l) {
		b[l+m] = (sg.y(l) *num) / den;
	
	  	if (l == 0) 
	    	sum_bN += b[m + l];
	  	else 
	    	sum_bN += 2.*b[m + l];

	  	if (l % r == 0)
	    	UppLog() << Format("\nb[%5d] = ", l) << FormatRational(b[l+M], 20);
	}
	 
	for (int l = 1; l < m+1; ++l) 
		b[m + l] =  b[m - l]; 
	
	UppLog() << "\nsumb = " << sum_bN.Simplify();
	VERIFY(sum_bN.Simplify() == 1.);
}

// val = 2/1 * 3/2 * 4/3 * ... If done n times, result has to be n
template<typename T>
T Loop() {
	T val = 1;
	for (T d = 1; d < T(100.); ++d) 
		val *= (d+1.)/d;
	return val;
}

void TestRational() {
	UppLog() << "\n\nRounding errors test";
	
	double dval = Loop<double>();
	Rational rval = Loop<Rational>();
	
	UppLog() << "\ndouble   == 100: " << ((dval == 100) ? "true" : "false");		// Fails
	UppLog() << "\nRational == 100: " << ((rval == 100.) ? "true" : "false");
	
	UppLog() << "\n";
	
	UppLog() << "\nsin() calculation";
	
	Polynomial<Rational> sinSeries;
	intInf fact = 1;
	int sign = 1;
	for (int i = 1; i < 25; i++) {
		fact *= i;
		if (!((i-1)%2)) {
			sinSeries[i] = Rational(intInf(sign), fact);
			sign = -sign;
		}
	}
	UppLog() << "\nsin() Taylor series is: " << sinSeries;
	
	Rational sin_1_3 = sinSeries.y(Rational(1, 3));	
	UppLog() << "\nsin(1/3) = " << sin_1_3;		
	UppLog() << "\nsin(1/3) = " << FormatRational(sin_1_3, 32);	
}
