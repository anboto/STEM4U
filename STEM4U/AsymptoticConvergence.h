// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2023, the Anboto author and contributors
#ifndef _STEM4U_AsymptoticConvergence_h_
#define _STEM4U_AsymptoticConvergence_h_

namespace Upp {


class AsymptoticConvergence {
public:
	double Get(const UVector<double> &y, double minr2) {
		ASSERT(x.size() == y.size());
		
		UVector<double> xx = clone(x);
		UVector<double> yy = clone(y);
		VectorXd ww = w;
		asy.SetWeight(w);
		
		while (xx.size() > 3) {
			asy.SetDegree(min(2, xx.size()-1)); 
		
			VectorXY data(xx, yy);
			double r2;
			asy.Fit(data, r2);
			if (!IsNull(r2) && r2 > minr2) 
				return asy.GetCoeff(0);

			// Removes the first element
			xx.Remove(0);		
			yy.Remove(0);
			VectorXd tmp = ww.tail(ww.size()-1);
			ww = tmp;
			double mx = ww.maxCoeff();
			ww /= mx;						// Normalized to the max value
			asy.SetWeight(ww);
		}
		return Null;
	}

	void SetX(const UVector<double> &_x) {			// Fitting weight is proportional to the area around each value
		x = clone(_x);
		w.resize(x.size());
		
		First(w) = x[1] - x[0];						// It is supposed the same area at the left of the first
		Last(w)  = x[x.size()-1] - x[x.size()-2];	// It is supposed the same area at the right of the last
		for (int i = 1; i < w.size()-1; ++i)
			w[i] = (x[i+1] - x[i-1])/2;				//	((w[i] - w[i-1]) + (w[i+1] - w[i]))/2
	
		double mx = w.maxCoeff();
		w /= mx;									// Normalized to the max value
	}
	AsymptoticEquation &GetEquation()			{return asy;}

private:
	AsymptoticEquation asy;
	UVector<double> x;
	VectorXd w;
};
 
}
 
#endif
