// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2023, the Anboto author and contributors
#ifndef _STEM4U_Wind_h_
#define _STEM4U_Wind_h_

double WindDensity(double wind, double rho, double rho0);

double RhoHeight(double height, double rho0, double height0, double temp0, double press0, double g = -1, double T0 = -1);
		
double GetWindShearCoeff(double height0, double height1, double wind0, double wind1, char shearLaw);

double GetWindHub(double height0, double heightHub, double wind0, char shearLaw, double coeff);

#endif
