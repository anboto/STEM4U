// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#ifndef _STEM4U_CrossCorrelation_h_
#define _STEM4U_CrossCorrelation_h_


namespace Upp {

template <class Range>
void XCorr(const Range &_x1, const Range &_x2, Range &r, Range &lags, char scale = 'n', size_t maxlag = 0) {
	ASSERT(_x2.size() == 0 || _x1.size() == _x2.size());
	
	size_t n = std::max(_x1.size(), _x2.size());
	ASSERT(maxlag <= n-1);

	if (maxlag == 0)
		maxlag = n-1;

	size_t M = size_t(PowInt(2., NextPow2(int(n + maxlag))));
	
	Eigen::VectorXd x1, x2;
	Copy(_x1, x1);
	
	if (_x2.size() == 0) {	
		Eigen::VectorXcd pre;
		Eigen::FFT<double> fft;
		
		PostPad(x1, M, 0.);
		fft.fwd(pre, x1);
		
		Eigen::VectorXcd post = pre.conjugate();
		
		for (int i = 0; i < pre.size(); ++i)
			pre[i] *= post[i];
		
		Eigen::VectorXd cor;
		fft.inv(cor, pre);
		
		Resize(r, 2*maxlag+1);
		r << cor.tail(maxlag), cor.head(maxlag+1);		
	} else {			
		Eigen::VectorXcd pre, post;
		Eigen::FFT<double> fft;
		
		M = max(M, x1.size() + maxlag);
		PrePad(x1, x1.size() + maxlag, 0.);
		PostPad(x1, M, 0.);
		fft.fwd(pre, x1);
		
		Copy(_x2, x2);
		
		PostPad(x2, M, 0.);
		fft.fwd(post, x2);
		
		post = post.conjugate();
		
		for (int i = 0; i < pre.size(); ++i)
			pre[i] *= post[i];
		
		Eigen::VectorXd cor;
		fft.inv(cor, pre);
		r = cor.segment(0, 2*maxlag + 1);
	}
	
	double dN = double(n);
	double dmaxlag = double(maxlag);
	
  	if (scale == 'b')				// biased	
    	r /= dN;
  	else if (scale == 'u') {		// unbiased
  		Vector<double> left, right;
  		LinSpaced(left, int(maxlag), dN-dmaxlag, dN-1);
  		left << dN;
  		LinSpaced(right, int(maxlag), dN-1, dN-dmaxlag);
  		left.Append(right);
    	r = r.cwiseQuotient(Eigen::Map<Eigen::VectorXd>(left, left.size()));
  	} else if (scale == 'c') {		// coeff
	    if (x2.size() == 0)
	      	r /= r[maxlag];
	    else {
	        double den = ::sqrt(x1.squaredNorm()*x2.squaredNorm());
	        if (den < 1e-10)
	            Zero(r);
	        else
	      		r /= den;
	    }
 	} else if (scale == 'n')		// none
 		;	
 	else 
 		NEVER_("Unknown scale");
	
    LinSpaced(lags, 2*int(maxlag)+1, -dmaxlag, dmaxlag);
}

template <class Range>
void XCorr(const Range &x, Range &r, Range &lags, char scale = 'n', size_t maxlag = 0) {
	Range y;
	XCorr(x, y, r, lags, scale, maxlag);	 
}

}

#endif
