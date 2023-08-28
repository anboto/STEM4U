// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2023, the Anboto author and contributors
#ifndef _STEM4U_AsymptoticConvergence_h_
#define _STEM4U_AsymptoticConvergence_h_

namespace Upp {


class AsymptoticConvergence {
public:
	double Get(const UVector<double> &y, double ratioY, double minr2, double rangeRatio, int way) {
		ASSERT(x.size() == y.size());
		
		xx = clone(x);
		yy = clone(y);
		VectorXd ww;	Copy(w, ww);
		
		// Remove outliers from avg with ratioY
		double avg = Avg(y);
		for (int i = yy.size()-1; i >= 0; --i) {
			if (abs(yy[i] - avg)/avg > ratioY ) {
				xx.Remove(i);		
				yy.Remove(i);
				if (ww.size() > 0)
					Remove(w, i);
			}
		}
		
		double rangex = Last(xx) - First(xx);
		double first = Last(xx) - rangex*rangeRatio;
		int idfirst;
		for (int i = xx.size()-1; i >= 0; --i) {
			if (xx[i] <= first) {
				idfirst = i;
				break;
			}
		}
		int num = max(4, xx.size() - idfirst);		
		
		AsymptoticEquation asy0;
		
		xx = Right(xx, num);
		yy = Right(yy, num);
		if (ww.size() > 0) {
			ww = Right(ww, num);
			asy0.SetWeight(ww);
		}
		
		VectorXY data(xx, yy);
		
		if (way == 0) {		// Get the degree that best fits
			double bestR2 = 0;
			for (int deg = 2; deg <= min(8, num); ++deg) {	
				asy0.SetDegree(deg); 
				
				double r2;
				asy0.Fit(data, r2);
				if (!IsNull(r2) && (bestR2 == 0 || r2 > bestR2)) {
					asy = asy0;
					bestR2 = r2;
				}
			}
			if (bestR2 != 0 && bestR2 >= minr2) 
				return asy.GetCoeff(0);
		} else {
			for (double fact = 0; fact > -0.2; fact -= 0.01) {
				UVector<double> xxx = clone(xx);
				UVector<double> yyy = clone(yy);
				VectorXd www;	Copy(ww, www);
				asy0.SetWeight(www);
			
				while (xxx.size() > 3) {						
					asy0.SetDegree(min(2, xxx.size()-1)); 	
				
					VectorXY data(xxx, yyy);
					double r2;
					asy0.Fit(data, r2);
					if (!IsNull(r2) && r2 > minr2 + fact) {
						asy = asy0;
						return asy.GetCoeff(0);
					}
					// Removes the first element
					xxx.Remove(0);		
					yyy.Remove(0);
					VectorXd tmp = www.tail(www.size()-1);
					www = tmp;
					double mx = www.maxCoeff();
					www /= mx;						// Normalized to the max value
					asy0.SetWeight(www);
				}
			}
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
	}
	AsymptoticEquation &GetEquation()			{return asy;}

	UVector<double> &GetRealX() 				{return xx;}
	UVector<double> &GetRealY() 				{return yy;}

private:
	AsymptoticEquation asy;
	UVector<double> x;
	UVector<double> xx, yy; 
	VectorXd w;
};
 
}
 
#endif
