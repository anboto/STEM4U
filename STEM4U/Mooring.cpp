// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#include <Core/Core.h>
#include <Eigen/Eigen.h>
#include "Sundials.h"
#include "Mooring.h"

#include <Functions4U/EnableWarnings.h>

namespace Upp {

MooringStatus Catenary(double rho_m, double rho_m3, double rho_water, double moorlen, double BL,
			  double xanchorvessel, double zanchor, double zvessel, 
			  double &Fhanchorvessel, double &Fvanchor, double &Fvvessel, 
			  double &xonfloor, Vector<double> &x, Vector<double> &z, int num) {
	if (zanchor < 0 || zvessel < 0)
		throw Exc("Catenary. Anchor and vessel z positions have to be positive");
	if (xanchorvessel < 0)
		throw Exc("Catenary. x distance between anchor and vessel has to be positive");
	
	const double g = 9.81;       		

	double lmb = rho_m*(rho_m3 - rho_water)/rho_m3;
	
	MooringStatus status;
	double straightLen = std::sqrt(sqr(xanchorvessel) + sqr(zanchor - zvessel));
	if (moorlen*1.02 < straightLen) {
		status = BROKEN;
		if (!IsNull(rho_m)) 
			Fhanchorvessel = Fvanchor = Fvvessel = xonfloor = 0;
		else
			Fhanchorvessel = Fvanchor = Fvvessel = xonfloor = Null;
	} else if (moorlen*1.02 > straightLen && moorlen < straightLen) {		// 2% elongation allowed
		status = TAUT;
		Fhanchorvessel = Fvanchor = Fvvessel = xonfloor = Null;				// Unknown without EA
		x.SetCount(num);
		z.SetCount(num);
		for (int i = 0; i < num; ++i) {
			x[i] = i*(xanchorvessel/(num-1.));
			z[i] = zanchor + i*((zvessel - zanchor)/(num-1.));
		}
    } else if (xanchorvessel < moorlen - zanchor - zvessel) {				
			status = LOOSE_ON_FLOOR;
			
			if (!IsNull(rho_m)) {
		        Fhanchorvessel = 0;
		        Fvanchor = lmb*g*zanchor;
		        Fvvessel = lmb*g*zvessel;
			} else
				Fhanchorvessel = Fvanchor = Fvvessel = Null;
			
	     	xonfloor = xanchorvessel;
	     	
			x.SetCount(4);
			z.SetCount(4);
			x[0] = 0;				z[0] = zanchor;
			x[1] = 0;				z[1] = 0;
			x[2] = xanchorvessel;	z[2] = 0;
			x[3] = xanchorvessel;	z[3] = zvessel;
    } else {
	  	Buffer<double> udata(1, 1);	
	  	bool dn = true;
	  	try {
			SolveNonLinearEquationsSun(udata, 1, [&](const double y[], double *residuals)->bool {
				double Blimit = y[0];
				double sqanchor = zanchor*(zanchor + 2*Blimit);
				if (sqanchor < 0)
					return false;
				double sqvessel = zvessel*(zvessel + 2*Blimit);
				if (sqvessel < 0)
					return false;
				residuals[0] = std::sqrt(sqanchor) + std::sqrt(sqvessel) - moorlen; 
				return true;
			});
	  	} catch (...) {
	  		dn = false;
	  	}
	  	if (!dn) {
	  		status = CALCULATION_PROBLEM;
			if (!IsNull(rho_m)) 
				Fhanchorvessel = Fvanchor = Fvvessel = xonfloor = 0;
			else
				Fhanchorvessel = Fvanchor = Fvvessel = xonfloor = Null;
	  	} else {
			double Blimit = udata[0];
				
		    double xanchorvessel_limit = Blimit*(acosh(zanchor/Blimit + 1) + acosh(zvessel/Blimit + 1));
			
			double delta = 0;
			if (num > 1)
				delta = xanchorvessel/(num - 1);
			x.SetCount(num);
			z.SetCount(num);
			
			if (xanchorvessel < xanchorvessel_limit) {
		        status = CATENARY_ON_FLOOR;
			
			  	Buffer<int> consdata(1, 2);		// B > 0
			  	
			  	Buffer<double> udata2(1, Blimit/2.);	
			  	bool done = true;
			  	try {
					SolveNonLinearEquationsSun(udata2, 1, [&](const double y[], double *residuals)->bool {
						double B = y[0];
						residuals[0] = xanchorvessel - B*acosh(zanchor/B + 1) - B*acosh(zvessel/B + 1) + std::sqrt(zanchor*(zanchor + 2*B)) + std::sqrt(zvessel*(zvessel + 2*B)) - moorlen; 
						return true;
					}, consdata);
			  	} catch(...) {
			  		done = false;
			  	}
			  	if (!done) {
			  		status = CALCULATION_PROBLEM;
					if (!IsNull(rho_m)) 
						Fhanchorvessel = Fvanchor = Fvvessel = xonfloor = 0;
					else
						Fhanchorvessel = Fvanchor = Fvvessel = xonfloor = Null;
			  	} else {
					double B = udata2[0];
			
			        double xcatanchor = B*acosh(zanchor/B + 1);
			        double xcatvessel = B*acosh(zvessel/B + 1);
					xonfloor = xanchorvessel - xcatanchor - xcatvessel;
			
					if (!IsNull(rho_m)) {
				        Fhanchorvessel = lmb*g*B;        
				        Fvanchor = Fhanchorvessel*sinh(xcatanchor/B); 
				        Fvvessel = Fhanchorvessel*sinh(xcatvessel/B); 
					} else
						Fhanchorvessel = Fvanchor = Fvvessel = Null;
			        
			        for (int i = 0; i < x.size(); ++i) {
			            x[i] = i*delta;
			    		if (x[i] < xcatanchor)
			        		z[i] = B*(cosh((xcatanchor - x[i])/B) - 1);
			            else if (x[i] > (xanchorvessel - xcatvessel))
			                z[i] = B*(cosh((max(0., x[i] - (xanchorvessel - xcatvessel)))/B) - 1);
			            else
			                z[i] = 0;
			        }
			  	}
		    } else if (xanchorvessel < std::sqrt(sqr(moorlen) - sqr(zvessel - zanchor))) {
		        status = CATENARY;
		        
				int neq = 2;
		        double deltaz = zvessel - zanchor;
		
			  	Buffer<double> udata2((size_t)neq);
			  	Buffer<int> consdata(2);
			  	consdata[0] = 2;			// B > 0
			  	consdata[1] = 0;
			  	
			  	double x1_0, B_0, x1, B;
			  	bool done = false;
			  	for (x1_0 = 0; x1_0 < abs(xanchorvessel) && !done; x1_0 += abs(xanchorvessel)/4) {
			  		for (B_0 = abs(Blimit); B_0 < 10*B_0 && !done; B_0 += 5*abs(Blimit)/4) {		
				  		udata2[0] = B_0;	
				  		udata2[1] = x1_0;
				  		bool don = true;
				  		try {
							SolveNonLinearEquationsSun(udata2, neq, [&](const double y[], double *residuals)->bool {
								double B = y[0];
								double x1 = y[1];
								
								residuals[0] = B*(sinh((x1 + xanchorvessel)/B) - sinh(x1/B)) - moorlen;
								residuals[1] = B*(cosh((x1 + xanchorvessel)/B) - cosh(x1/B)) - deltaz;
								
								return true;
							}, consdata);
				  		} catch (...) {
				  			don = false;
				  		}
				  		if (don) {
							B = udata2[0];
							x1 = udata2[1];		
							if (abs((x1 + xanchorvessel)/B) < 3)
								done = true;
				  		}
			  		}
			  	}
			  	if (!done) {
					//throw Exc("Catenary. Solving problem");
					status = CALCULATION_PROBLEM;
					if (!IsNull(rho_m)) 
						Fhanchorvessel = Fvanchor = Fvvessel = xonfloor = 0;
					else
						Fhanchorvessel = Fvanchor = Fvvessel = xonfloor = Null;
			  	} else {
					if (!IsNull(rho_m)) {
				        Fhanchorvessel = lmb*g*abs(B);    
				        Fvanchor = -Fhanchorvessel*sinh(x1/B);  			
				        Fvvessel = Fhanchorvessel*sinh((x1 + xanchorvessel)/B);  	
					} else
						Fhanchorvessel = Fvanchor = Fvvessel = Null;
					
			        xonfloor = 0;
			        
			        for (int i = 0; i < x.size(); ++i) {
			            x[i] = i*delta;
			        	z[i] = B*(cosh((x[i] + x1)/B) - 1) - B*(cosh(x1/B) - 1) + zanchor;
			        }
			  	}
		    } else {
				status = CALCULATION_PROBLEM;
				if (!IsNull(rho_m)) 
					Fhanchorvessel = Fvanchor = Fvvessel = xonfloor = 0;
				else
					Fhanchorvessel = Fvanchor = Fvvessel = xonfloor = Null;
		    }
    	}
    }
    
	if (!IsNull(BL))    
    	if (max(std::sqrt(sqr(Fhanchorvessel) + sqr(Fvanchor)), std::sqrt(sqr(Fhanchorvessel) + sqr(Fvvessel))) >= BL)
        	status = BL_EXCEDEED;
    
    if (status == BL_EXCEDEED || status == BROKEN || status == CALCULATION_PROBLEM) {
        x.SetCount(2);
        z.SetCount(2);
        x[0] = 0;				z[0] = zanchor;
		x[1] = xanchorvessel;	z[1] = zvessel;
    }
        
	return status;
}

MooringStatus Catenary(double rho_m, double rho_m3, double rho_water, double moorlen, double BL,
			  double xanchorvessel, double zanchor, double zvessel, 
			  double &Fhanchorvessel, double &Fvanchor, double &Fvvessel, double &xonfloor) {
	Vector<double> x, z;
				      
	return Catenary(rho_m, rho_m3, rho_water, moorlen, BL, xanchorvessel, 
			zanchor, zvessel, Fhanchorvessel, Fvanchor, Fvvessel, xonfloor, x, z, 0);				      
}

MooringStatus Catenary(double moorlen, double xanchorvessel, double zanchor, double zvessel, double &xonfloor,
				Vector<double> &x, Vector<double> &z, int num) {
	double rho_m = Null,  rho_m3 = Null, rho_water = Null, BL = Null, Fhanchorvessel, Fvanchor, Fvvessel;		      

	return Catenary(rho_m, rho_m3, rho_water, moorlen, BL, xanchorvessel, 
			zanchor, zvessel, Fhanchorvessel, Fvanchor, Fvvessel, xonfloor, x, z, num);
}

MooringStatus Catenary(double moorlen, double xanchorvessel, double zanchor, double zvessel, double &xonfloor) {
	double rho_m = Null,  rho_m3 = Null, rho_water = Null, BL = Null, Fhanchorvessel, Fvanchor, Fvvessel;		      
	Vector<double> x, z;
	return Catenary(rho_m, rho_m3, rho_water, moorlen, BL, xanchorvessel, 
			zanchor, zvessel, Fhanchorvessel, Fvanchor, Fvvessel, xonfloor, x, z, 0);
}


bool CatenaryGetLen0(double xonfloor, double xanchorvessel, double zanchor, double zvessel, 
					double &moorlen) {
	if (xonfloor >= xanchorvessel) {
		moorlen = xanchorvessel + 2*(xonfloor - xanchorvessel) + zanchor + zvessel;
		return false;
	} 
	
	Buffer<int> consdata(2, 2);		// B > 0 && moorlen > 0
  	Buffer<double> udata(2);	
  	udata[0] = 1;
  	udata[1] = std::sqrt(sqr(xanchorvessel) + sqr(zanchor - zvessel));
	SolveNonLinearEquationsSun(udata, 2, [&](const double y[], double *residuals)->bool {
		double B = y[0];
		double moorlen = y[1];
		residuals[0] = xanchorvessel - B*acosh(zanchor/B + 1) - B*acosh(zvessel/B + 1) + std::sqrt(zanchor*(zanchor + 2*B)) + std::sqrt(zvessel*(zvessel + 2*B)) - moorlen; 
    	double xcatanchor = B*acosh(zanchor/B + 1);
    	double xcatvessel = B*acosh(zvessel/B + 1);
		residuals[1] = xanchorvessel - xcatanchor - xcatvessel - xonfloor;
		return true;
	}, consdata);		
	//double B = udata[0];
	moorlen = udata[1];

	return true;
					}

MooringStatus CatenaryGetLen(double xonfloor, double xanchorvessel, double zanchor, double zvessel, 
			double &moorlen) {
	if (!CatenaryGetLen0(xonfloor, xanchorvessel, zanchor, zvessel, moorlen))
		return LOOSE_ON_FLOOR;
	
	return Catenary(moorlen, xanchorvessel, zanchor, zvessel, xonfloor);
}

MooringStatus CatenaryGetLen(double rho_m, double rho_m3, double rho_water, double xonfloor, 
			double BL, double xanchorvessel, double zanchor, double zvessel, 
			double &Fhanchorvessel, double &Fvanchor, double &Fvvessel, double &moorlen) {
	if (!CatenaryGetLen0(xonfloor, xanchorvessel, zanchor, zvessel, moorlen))
		return LOOSE_ON_FLOOR;
	
	return Catenary(rho_m, rho_m3, rho_water, moorlen, BL, xanchorvessel, zanchor, zvessel, 
					Fhanchorvessel, Fvanchor, Fvvessel, xonfloor);
}

const char *MooringStatusStr(MooringStatus status) {
	const char *str[7] = {t_("loose on seabed"), t_("catenary on seabed"), t_("catenary"), t_("taut"),
						  t_("line length exceeded"), t_("break load exdeeded"), t_("Calculation problem")};
	ASSERT(int(status) < 6);
	return str[int(status)];
}

bool IsOK(MooringStatus status) { 
	return status < BROKEN;
}

}