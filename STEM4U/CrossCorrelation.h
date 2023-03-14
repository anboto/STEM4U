// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#ifndef _STEM4U_CrossCorrelation_h_
#define _STEM4U_CrossCorrelation_h_


namespace Upp {

template <class Range>
void XCorr(const Range &_x, const Range &_y, Range &r, Range &lags, char scale = 'n', size_t maxlag = 0) {
	ASSERT(_y.size() == 0 || _x.size() == _y.size());
	
	size_t n = std::max(_x.size(), _y.size());

	if (maxlag == 0)
		maxlag = n-1;
	
	ASSERT(maxlag <= n-1);

	size_t P = _x.size();
	size_t M = size_t(pow(2, NextPow2(int(n + maxlag))));
	
	Eigen::VectorXd x, y;
	
	Copy(_x, x);
	
	if (_y.size() == 0) {	
		Eigen::VectorXcd pre;
		Eigen::FFT<double> fft;
		
		PostPad(x, M, 0.);
		fft.fwd(pre, x);
		
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
		
		PrePad(x, x.size() + maxlag, 0.);
		PostPad(x, M, 0.);
		fft.fwd(pre, x);
		
		Copy(_y, y);
		
		PostPad(y, M, 0.);
		fft.fwd(post, y);
		
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
	    if (y.size() == 0)
	      	r /= r[maxlag];
	    else {
	        double den = ::sqrt(x.squaredNorm()*y.squaredNorm());
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
