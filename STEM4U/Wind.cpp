// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2023, the Anboto author and contributors
#include <Core/Core.h>

#include "Wind.h"

#include <Functions4U/EnableWarnings.h>

namespace Upp {

// Svenningsen correction (WindPro y OpenWind) following ITC
double WindDensity(double wind, double rho, double rho0) {
    double gamma;
    
    if (wind <= 8)
        gamma = 1 / 3.;
    else if (wind >= 12)
        gamma = 2 / 3.;
    else
        gamma = 1 / 3. + (wind - 8) * (2 / 3. - 1 / 3.) / (12 - 8);

    return wind * pow(rho / rho0, gamma);
}

// https://en.wikipedia.org/wiki/Density_of_air & https://www.homerenergy.com/products/pro/docs/latest/altitude.html
double RhoHeight(double height, double rho0, double height0, double temp0, double press0, double g, double T0) {
	const double R = 8.3144598;		// J/mol/K	Ideal gas constant
	const double M = 0.0289654; 	// kg/mol	Molar mass of dry air
	const double L = 0.0065;		// K/m		Temperature lapse rate
	
	if (g < 0)
		g = 9.80665;		// m/s2
	if (T0 < 0)
		T0= 273.15;	

	if (rho0 < 0) {
		if (temp0 < 0 || press0 < 0) 
			rho0 = 1.225;
		else
			rho0 = press0*M/R/(temp0 + T0);
	}
	if (temp0 < 0 || press0 < 0) 
		return rho0*exp(-(height - height0)/10400);	
	else {
		double temp = temp0 - L*(height - height0);
		double press = press0*pow((1 - L*height/(temp0 + T0)), g*M/R/L);
		return rho0*press/press0*((temp0 + T0)/(temp + T0));
	}	
}
		
double GetWindShearCoeff(double height0, double height1, double wind0, double wind1, char shearLaw) {
	if (wind0 == 0 || wind1 == 0 || height1 == 0)
		return -1;
	
	if (shearLaw == 'e') 			// Exponential, coeff. alpha
		return log(wind0/wind1)/log(height0/height1);
	else 						// Logarithmic, coeff. z0
		return exp((wind1*log(height0) - wind0*log(height1))/(wind1 - wind0));
}

double GetWindHub(double height0, double heightHub, double wind0, int shearLaw, double coeff) {
	if (shearLaw == 'e') 			// Exponential, coeff. alpha
		return wind0*pow(heightHub/height0, coeff);
	else 						// Logarithmic, coeff. z0
		return wind0*(1 + log(heightHub/height0)/log(height0/coeff));
}

}