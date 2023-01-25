// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#ifndef _STEM4U_CrossCorrelation_h_
#define _STEM4U_CrossCorrelation_h_


namespace Upp {

template <class Range>
void XCorr(const Range &_x, const Range &_y, Range &_r, Range &lags, char scale = 'n', size_t maxlag = 0) {
	size_t N = std::max(_x.size(), _y.size());

	if (maxlag == 0)
		maxlag = N-1;
	
	ASSERT(maxlag <= N-1);

	size_t P = _x.size();
	size_t M = size_t(pow(2, NextPow2(int(N + maxlag))));
	
	Eigen::VectorXd x, y, r;
	
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
	
	double dN = double(N);
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
	    else
	      	r /= ::sqrt(x.squaredNorm()*y.squaredNorm());
 	} else if (scale == 'n')		// none
 		;	
 	else 
 		NEVER_("Unknown scale");
	
	Copy(r, _r);
	
    LinSpaced(lags, 2*int(maxlag)+1, -dmaxlag, dmaxlag);
}

template <class Range>
void XCorr(const Range &x, Range &r, Range &lags, char scale = 'n', size_t maxlag = 0) {
	Range y;
	XCorr(x, y, r, lags, scale, maxlag);	 
}

}

#endif
