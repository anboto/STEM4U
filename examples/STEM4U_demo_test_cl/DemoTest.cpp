// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#include <Core/Core.h>
#include <Functions4U/Functions4U.h>
#include <Eigen/Eigen.h>
#include <STEM4U/IntInf.h>
#include <STEM4U/Polynomial.h>
#include <STEM4U/Sundials.h>
#include <STEM4U/Integral.h>
#include <STEM4U/TravellingSalesman.h>
#include <STEM4U/ShortestPath.h>
#include <STEM4U/SeaWaves.h>
#include <STEM4U/Utility.h>
#include <STEM4U/CrossCorrelation.h>
#include <ScatterDraw/DataSource.h>
#include <STEM4U/Rootfinding.h>
#include <Eigen/MultiDimMatrix.h>
#include <Painter/Painter.h>
#include <plugin/png/png.h>

using namespace Upp;
using namespace Eigen;


void TestVectorMatrixHealing() {
	UppLog() << "\n\nUncomplete vector and matrix healing and filling";
	
	{
		VectorXd x(4), nx;
		VectorXd z(4), nz;
		
		x << 2, 4, 6, 10;
		z << 0, 1, 2, 4;
		
		UppLog() << "\nInitial data";
		UppLog() << "\nx:\n" << x.transpose();
		UppLog() << "\nz:\n" << z.transpose();
		
		GapFilling(x, z, nx, nz, true, 10);
		
		UppLog() << "\nFilled with zero";
		UppLog() << "\nx:\n" << nx.transpose();
		UppLog() << "\nz:\n" << nz.transpose();
	
		GapFilling(x, z, nx, nz, false, 10);
		
		UppLog() << "\nInterpolated";
		UppLog() << "\nx:\n" << nx.transpose();
		UppLog() << "\nz:\n" << nz.transpose();
		
		UppLog() << "\n";
	}
	{
		VectorXd x(4), y(3), nx, ny;
		MatrixXd z(3, 4), nz;
		
		x << 2, 4, 6, 10;
		y << 1, 3, 7;
		z << 0, 1, 2, 4,
			 1, 2, 3, 5,
			 3, 4, 5, 7;
		
		UppLog() << "\nInitial data";
		UppLog() << "\nx:\n" << x.transpose();
		UppLog() << "\ny:\n" << y.transpose();
		UppLog() << "\nz:\n" << z;
		
		GapFilling(x, y, z, nx, ny, nz, true, 10);
		
		UppLog() << "\nFilled with zero";
		UppLog() << "\nx:\n" << nx.transpose();
		UppLog() << "\ny:\n" << ny.transpose();
		UppLog() << "\nz:\n" << nz;
	
		GapFilling(x, y, z, nx, ny, nz, false, 10);
		
		UppLog() << "\nInterpolated";
		UppLog() << "\nx:\n" << nx.transpose();
		UppLog() << "\ny:\n" << ny.transpose();
		UppLog() << "\nz:\n" << nz;
		
		UppLog() << "\n";
	}
}

void TestRootfinding() {
	UppLog() << "\n\nRoot finding functions";
	
	UppLog() << "\n'Root polishing':";
	{
		int nf = 0, ndf = 0;
		double y = NewtonRaphson<double>([&](double x)->double {nf++; return sin(x);}, [&](double x)->double {ndf++; return cos(x);}, 3*M_PI/4, 0.001, 50);
		UppLog() << Format("\nNewtonRaphson:      %2d f(x) %2d df(x) y = %f", nf, ndf, y); 
		VERIFY(abs(y - M_PI) < 0.001);
	}
	{
		int nf = 0;
		double y = QuasiNewtonRaphson<double>([&](double x)->double {nf++; return sin(x);}, 3*M_PI/4, 0.01, 0.001, 50);
		UppLog() << Format("\nQuasiNewtonRaphson: %2d f(x)          y = %f", nf, y);
		VERIFY(abs(y - M_PI) < 0.001);
	}
	
	UppLog() << "\n'Root bracketing':";
	
	double from, to;
	RootBracketing<double>([&](double x)->double {return sin(x);}, 3*M_PI/4, 0, 2*M_PI, M_PI/4, from, to, false);
	{
		int nf = 0;
		double y = Bisection<double>([&](double x)->double {nf++; return sin(x);}, from, to, 0.001, 50);
		UppLog() << Format("\nBisection:          %2d f(x)          y = %f", nf, y);
		VERIFY(abs(y - M_PI) < 0.001);
	}
	{
		int nf = 0;
		double y = Brent<double>([&](double x)->double {nf++; return sin(x);}, from, to, 0.001, 50);
		UppLog() << Format("\nBrent:              %2d f(x)          y = %f", nf, y);
		VERIFY(abs(y - M_PI) < 0.001);
	}
		
	UppLog() << "\n";
}

void TestIntInf() {
	UppLog() << "\n\nintInf demo. A signed integer type with arbitrary-precision including the usual arithmetic.";
	
	intInf a = "12345678901234567890";
	intInf b = 2, c;
	
	c = a%b;		UppLog() << "\na%%b:      "	<< c;
	
	VERIFY(c == 0);
	
	c = a;
	c += b;			UppLog() << "\nc += b:    " << c;
	c -= b;			UppLog() << "\nc -= b:    " << c;
	c = c + 2;		UppLog() << "\nc = 2 + b: " << c;
	c = c - 2;		UppLog() << "\nc = 2 - b: " << c;
	c *= 2;			UppLog() << "\nc *= b:    " << c;
	c /= 2;			UppLog() << "\nc /= b:    " << c;
	c = c * 2;		UppLog() << "\nc = c * 2: " << c;
	c = c / 2;		UppLog() << "\nc = c / 2: " << c;

	VERIFY(c == a);
}

void TestTravellingSalesman() {
	UppLog() << "\n\nTravelling salesman";

	const UVector<Point_<int>> points = {{0, 0},{4, 4},{4, 0},{2, 4},{0, 4},{4, 2},{0, 2},{2, 0}};

	UVector<int> orderp;
	int distp = TSP(points, orderp, TSP_NEAREST_NEIGHBOR);
	UppLog() << "\nTotal distance between points is: " << distp;
	VERIFY(distp == 16);
	String sorderp;
	for (int i = 0; i < orderp.size(); ++i) {
		if (i > 0)
			sorderp << " -> ";
		sorderp << orderp[i];
	}
	UppLog() << "\nOrder is: " << sorderp;
	UppLog() << "\n";
	VERIFY(sorderp == "0 -> 7 -> 2 -> 5 -> 1 -> 3 -> 4 -> 6 -> 0");
		
	// Example from https://developers.google.com/optimization/routing/tsp#printer
	const UVector<UVector<int>> cities = {
		{0, 2451, 713, 1018, 1631, 1374, 2408, 213, 2571, 875, 1420, 2145, 1972},
		{2451, 0, 1745, 1524, 831, 1240, 959, 2596, 403, 1589, 1374, 357, 579},
		{7133, 1745, 0, 355, 920, 803, 1737, 851, 1858, 262, 940, 1453, 1260},
		{1018, 1524, 355, 0, 700, 862, 1395, 1123, 1584, 466, 1056, 1280, 987},
		{1631, 831, 920, 700, 0, 663, 1021, 1769, 949, 796, 879, 586, 371},
		{1374, 1240, 803, 862, 663, 0, 1681, 1551, 1765, 547, 225, 887, 999},
		{2408, 959, 1737, 1395, 1021, 1681, 0, 2493, 678, 1724, 1891, 1114, 701},
		{213, 2596, 851, 1123, 1769, 1551, 2493, 0, 2699, 1038, 1605, 2300, 2099},
		{2571, 403, 1858, 1584, 949, 1765, 678, 2699, 0, 1744, 1645, 653, 600},
		{875, 1589, 262, 466, 796, 547, 1724, 1038, 1744, 0, 679, 1272, 1162},
		{1420, 1374, 940, 1056, 879, 225, 1891, 1605, 1645, 679, 0, 1017, 1200},
		{2145, 357, 1453, 1280, 586, 887, 1114, 2300, 653, 1272, 1017, 0, 504},
		{1972, 579, 1260, 987, 371, 999, 701, 2099, 600, 1162, 1200, 504, 0},
	};

	UVector<int> order;
	int dist = TSP(cities, order, TSP_CONSECUTIVE);
	
	UppLog() << "\nTotal distance between cities is: " << dist;
	VERIFY(dist == 7293);
	String sorder;
	for (int i = 0; i < order.size(); ++i) {
		if (i > 0)
			sorder << " -> ";
		sorder << order[i];
	}
	UppLog() << "\nOrder is: " << sorder;
	VERIFY(sorder == "0 -> 7 -> 2 -> 3 -> 4 -> 12 -> 6 -> 8 -> 1 -> 11 -> 10 -> 5 -> 9 -> 0");
}

void TestShortestPath() {
	UppLog() << "\n\nShortest path";
	
	using Seg = SegSP<int>;
	
	UVector<UVector<Seg>> adjList;
	
	adjList.Add() << Seg(1, 2) << Seg(2, 3);
    adjList.Add() << Seg(0, 2) << Seg(5, 1);
    adjList.Add() << Seg(0, 3) << Seg(5, 2);
    adjList.Add() << Seg(1, 4) << Seg(4, 1) << Seg(6, 2);
    adjList.Add() << Seg(3, 1) << Seg(5, 2) << Seg(6, 1);
    adjList.Add() << Seg(1, 1) << Seg(2, 2) << Seg(4, 2) << Seg(6, 2);
    adjList.Add() << Seg(3, 2) << Seg(4, 1) << Seg(5, 2);

	int start = 0, end = 6;
	
	UppLog() << "\n\nDijkstra method";
	{
		UVector<Seg> dist = Dijkstra(adjList, start);
		
		UppLog() << Format("\nShortest distance from node %d", start);
    	for(int i = 0; i < dist.size(); i++) 
    		UppLog() << Format("\nto node %d: %d", i, dist[i].weight);
		VERIFY(dist[end].weight == 5);
	    int currnode = end;
	    String sorder;
	    sorder << end;
	    while(currnode != start) {
	        currnode = dist[currnode].node;
	        sorder << " <- " << currnode;
	    }
		UppLog() << Format("\nThe path from %d to %d is: %s", start, end, sorder);
		VERIFY(sorder == "6 <- 5 <- 1 <- 0"); 
	}
	UppLog() << "\n\nBellman-Ford method";
	{
		UVector<Seg> dist = BellmanFord(adjList, start);
		
		UppLog() << Format("\nShortest distance from node %d", start);
    	for(int i = 0; i < dist.size(); i++) 
    		UppLog() << Format("\nto node %d: %d", i, dist[i].weight);
		VERIFY(dist[end].weight == 5);
	    int currnode = end;
	    String sorder;
	    sorder << end;
	    while(currnode != start) {
	        currnode = dist[currnode].node;
	        sorder << " <- " << currnode;
	    }
		UppLog() << Format("\nThe path from %d to %d is: %s", start, end, sorder);
		VERIFY(sorder == "6 <- 5 <- 1 <- 0"); 
	}
	UppLog() << "\n\nFloyd-Warshall method";
	{
		const int inf = std::numeric_limits<int>::max();
		UVector<UVector<int>> adjMatrix = {
			{0,  2,  3,  inf,inf,inf,inf},
			{2,  0,  inf,4,  inf,1,  inf},
			{3,  inf,0,  inf,inf,2,  inf},
			{inf,4,  inf,0,  1,  inf,2}, 
			{inf,inf,inf,1,  0,  2,  1},
			{inf,1  ,2,  inf,2,  0,  2},
			{inf,inf,inf,2,  1,  2,  0}
		};
		
		UVector<UVector<int>> dist = FloydWarshall(adjMatrix);

    	UppLog() << "\nShortest distance for all nodes";
    	for(int i = 0; i < dist.size(); i++)
        	for(int j = 0; j < dist.size(); j++)
            	UppLog() << "\nNode " << i << " to node " << j << ": " << dist[i][j];
    
	}
}

void TestDAESolver() {
	UppLog() << "\n\nSolveDAE() solves an harmonic oscillator m·d2x + k·x = 0"; 
	double y[]  = {2, 0};
	double dy[] = {0, 0};
	double m = 1, k = 0.5;
	SolveDAE(y, dy, 2, 0.1, 10, 
		[&](double t, Eigen::Index iiter, const double y[], const double dy[], double residual[])->int {
			residual[0] = m*dy[1] + k*y[0];
			residual[1] = y[1] - dy[0];
			return true;
		}, 2,
		[&](double t, Eigen::Index iiter, const double y[], const double dy[], double residual[])->int {
			residual[0] = y[0] - 0.0001;
  			residual[1] = y[1] - 0.0001;
			return true;
		},
		[&](double t, Eigen::Index iiter, const double y[], const double dy[], bool isZero, int *whichZero)->bool {
			UppLog() << Format("\n>T: %7.4f %8.4f %8.4f %s", t, y[0], y[1], isZero ? "Y" : "");
			return true;
		}
	);
}


void TestIntegral() {	
	UppLog() << "\n\nIntegral demo";
	
	double R = 1;
	double res = M_PI*sqr(R)/2;
	UppLog() << "\nSemicircle integral value is " << res;
	UppLog() << "\nNumerically integrated with simple and composite versions:";
	UppLog() << Format("\n%s\t%s\t%s\t%s\t%s\t\t%s\t%s", "#Points", "Trapezoidal", "Simpson 1/3", "Simpson 3/8", "Spline", "Hermite 3 point", "Hermite 5 point");
	for (int nump = 5; nump <= 30; ++nump) {
		double dx = 2*R/(nump-1);
		VectorXd y(nump), x(nump);
		for (int i = 0; i < nump; ++i) {
			x[i] = 2*R*i/(nump-1) - R;
			y[i] = ::sqrt(sqr(R) - sqr(x[i]));
		}
		
		double yx_trap    = Integral(x, y, TRAPEZOIDAL);
		double yx_simp13  = Integral(x, y, SIMPSON_1_3);
		double yx_simp38  = Integral(x, y, SIMPSON_3_8);
		double yx_spline  = Integral(x, y, SPLINE);
		double ydx_trap   = Integral(y, dx, TRAPEZOIDAL);
		double ydx_simp13 = Integral(y, dx, SIMPSON_1_3);
		double ydx_simp38 = Integral(y, dx, SIMPSON_3_8);
		double ydx_herm3  = Integral(y, dx, HERMITE_3);
		double ydx_herm5  = Integral(y, dx, HERMITE_5);
		double ydx_spline = Integral(y, dx, SPLINE);
		
		UppLog() << Format("\n%d", nump);
		UppLog() << Format("\t%7.5f=%7.5f", yx_trap, ydx_trap);
		UppLog() << Format("\t%7.5f=%7.5f", yx_simp13, ydx_simp13);
		UppLog() << Format("\t%7.5f=%7.5f", yx_simp38, ydx_simp38);
		UppLog() << Format("\t%7.5f=%7.5f", yx_spline, ydx_spline);
		UppLog() << Format("\t%7.5f", ydx_herm3);
		UppLog() << Format("\t\t%7.5f", ydx_herm5);
		VERIFY(yx_trap - ydx_trap < 1E-10);			VERIFY(abs(yx_trap   - res)/res < 0.15);
		VERIFY(yx_simp13 - ydx_simp13 < 1E-10);		VERIFY(abs(yx_simp13 - res)/res < 0.15);
		VERIFY(yx_simp38 - ydx_simp38 < 1E-10);		VERIFY(abs(yx_simp38 - res)/res < 0.15);
		VERIFY(yx_spline - ydx_spline < 1E-10);		VERIFY(abs(yx_spline - res)/res < 0.15);
	}
}

void TestIntegralSinCos() {
	UppLog() << "\n\nIntegral sin cos demo";
	
	double a = M_PI;
	double b = 2*M_PI;
	double t = 8;
	double exact;
	
	exact = 0.0051257792400619;
	UppLog() << "\nintegral pi->2*pi 4*exp(-x/2)*cos(8*x) dx = " << exact << "\n#\tError";
	
	for (int it = 0; it < 10; it++) {
	    int n = 10 + it*20 + 1;
		double dx = (b-a)/(n-1);
	
	    VectorXd x(n), f(n);
	    for (int i = 0; i < n; i++ ) {
	        x[i] = a + i*dx;
	        f[i] = 4*exp(-x[i]/2);
	    }
	    
	    double res = IntegralSinCos(a, (b-a)/(n-1), f, t, true);
		double error = (res - exact)/exact;
	
	    UppLog() << "\n" << n << "\t" << error;
	    
	    if (n == 111) 
	        VERIFY(abs(error) < 1E-6);
	}
	
	exact = 0.082012467840991;
	UppLog() << "\nintegral pi->2*pi 4*exp(-x/2)*sin(8*x) dx = " << exact << "\n#\tError";
	
	for (int it = 0; it < 10; it++) {
	    int n = 10 + it*20 + 1;
		double dx = (b-a)/(n-1);
	
	    VectorXd x(n), f(n);
	    for (int i = 0; i < n; i++ ) {
	        x[i] = a + i*dx;
	        f[i] = 4*exp(-x[i]/2);
	    }
	    
	    double res = IntegralSinCos(a, (b-a)/(n-1), f, t, false);
		double error = (res - exact)/exact;
	
	    UppLog() << "\n" << n << "\t" << error;
	    
	    if (n == 111) 
	        VERIFY(abs(error) < 1E-9);
	}
}
	
void TestSeaWaves() {
	UppLog() << "\n\nSeaWaves demo";

	double T = 12;
	double depth = 50;
	double H = 2;
	
	double rho = 1025, g = 9.81;
	
	double waveNumber  = SeaWaves::WaveNumber(T, depth, g, false);
	UppLog() << "\n" << Format("Wave number: %f rad/m", waveNumber);
	double waveNumberE = SeaWaves::WaveNumber(T, depth, g, true);
	UppLog() << "\n" << Format("Wave number (exact): %f rad/m", waveNumberE);
	double waveLength = SeaWaves::WaveLength(T, depth, g);
	UppLog() << "\n" << Format("Wave length: %f m", waveLength);
	double c = SeaWaves::Celerity(T, depth, g);		
	UppLog() << "\n" << Format("Celerity: %f m/s", c);
	double gc = SeaWaves::GroupCelerity(T, depth, g);
	UppLog() << "\n" << Format("Group celerity: %f m/s", gc);
	SeaWaves::SEA_TYPE seaType = SeaWaves::GetSeaType(T, depth, g);
	UppLog() << "\n" << Format("Sea: %s", seaType == SeaWaves::SHALLOW ? "shallow" : seaType == SeaWaves::INTERMEDIATE ? "intermediate" : "deep");
	double power = SeaWaves::Power(T, H, depth, g, rho);
	UppLog() << "\n" << Format("Power: %f kW/m", power);
		
	double Tz = 12;
	double gamma = 2;
	double Tp = Tp_fTz(Tz, gamma);
	double gamma2 = gamma_fTp_Tz(Tp, Tz);
	UppLog() << "\n" << Format("Tp: %.2f, Tz: %.2f, gamma: %.4f, %.4f", Tp, Tz, gamma, gamma2);
	VERIFY(abs(gamma - gamma2) < 0.000001);

	double freq = 2*M_PI/T;
	double psd = SeaWaves::JONSWAP_Spectrum(H, T, 3.3, freq);
	UppLog() << "\n" << Format("JONSWAP PSD (%f): %f", freq, psd);

	{
		SeaWaves waves;
		
		double Hs = 2, Tp =	12, h = 70;
		waves.Init(Tp, Hs, 0, h);
		double x = 100, y = 100, z = -10, t = 10;
		waves.Calc(x, y, z, t);
		UppLog() << "\n" << Format("Sea data for Hs: %.2f m, Tp; %.2f s, at x: %.2f m, y: %.2f m, z: %.2f m, t: %.3f s", Hs, Tp, x, y, z, t);
		UppLog() << "\n" << Format("Free surface z: %f m = %f m", waves.zSurf, waves.ZSurf(x, y, t));
		VERIFY(abs(waves.zSurf - waves.ZSurf(x, y, t)) < 0.000001);
		VERIFY(abs(waves.zSurf + 0.0758260574) < 0.000001);
		UppLog() << "\n" << Format("vx: %f m/s, vy: %f m/s, vz: %f m/s", waves.vx, waves.vy, waves.vz);
		VERIFY(abs(waves.vz + 0.168306635) < 0.000001);
		UppLog() << "\n" << Format("ax: %f m/s2, ay: %f m/s2, az: %f m/s2", waves.ax, waves.ay, waves.az);
		VERIFY(abs(waves.az - 0.00269468) < 0.000001);
		UppLog() << "\n" << Format("p: %.3f Pa = %.3f Pa", waves.p, waves.Pressure(x, y, z, t));
		VERIFY(abs(waves.p - waves.Pressure(x, y, z, t)) < 0.000001);
		
		double duration = 60*60, deltaT = 0.1;
		VectorXd t_, et;
		VERIFY(waves.GetZSurf(t_, et, 0, 0, duration, deltaT));
		
		WaveParam param;
		SeaWaves::GetWaveParam(param, et, deltaT, h, g, rho);
		Cout() << Format("\nWave obtained %s", param.ToString());
	}
}


void TestXCorr() {
	VectorXd n, y;
	LinSpaced(n, 16, 0., 15.);
	VectorXd x = Eigen::pow(0.84, n.array());
	CircShift(x, 5, y); 
	
	VectorXd R, lags;
	
	XCorr(x, y, R, lags);
	UVector<double> realR = {0.1749,0.35513,0.54619,0.75389,0.98456,1.2452,1.5439,1.8896,2.2928,2.766,3.3234,2.8648,2.4935,2.1982,1.9699,1.8017,1.5026,1.2494,1.0343,0.85067,0.69298,0.55641,0.43679,0.33049,0.23425,0.14516,0.060494,0.046321,0.033559,0.02182,0.010746};
	UVector<double> realLags = {-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
	VERIFY(CompareDecimals(R, realR, 2));
	VERIFY(CompareDecimals(lags, realLags, 4));

	XCorr(x, y, R, lags, 'b');
	realR = {0.010931,0.022196,0.034137,0.047118,0.061535,0.077828,0.096492,0.1181,0.1433,0.17287,0.20771,0.17905,0.15585,0.13739,0.12312,0.1126,0.093915,0.078089,0.064643,0.053167,0.043311,0.034775,0.027299,0.020655,0.014641,0.0090727,0.0037809,0.002895,0.0020974,0.0013638,0.00067165};
	VERIFY(CompareDecimals(R, realR, 3));
	VERIFY(CompareDecimals(lags, realLags, 5));

	XCorr(x, y, R, lags, 'u');
	realR = {0.1749,0.17757,0.18206,0.18847,0.19691,0.20754,0.22055,0.2362,0.25476,0.2766,0.30213,0.23874,0.19181,0.15702,0.13133,0.1126,0.10018,0.089245,0.079561,0.070889,0.062998,0.055641,0.048532,0.041311,0.033465,0.024194,0.012099,0.01158,0.011186,0.01091,0.010746};
	VERIFY(CompareDecimals(R, realR, 3));
	VERIFY(CompareDecimals(lags, realLags, 4));

	XCorr(x, y, R, lags, 'c');
	realR = {0.051686,0.10495,0.16141,0.22279,0.29095,0.36799,0.45624,0.55839,0.67757,0.81739,0.98212,0.8466,0.73688,0.64961,0.58214,0.53242,0.44405,0.36922,0.30565,0.25139,0.20479,0.16443,0.12908,0.097664,0.069226,0.042898,0.017877,0.013688,0.0099172,0.0064482,0.0031757};
	VERIFY(CompareDecimals(R, realR, 3));
	VERIFY(CompareDecimals(lags, realLags, 4));

	XCorr(x, R, lags);
	realR = {0.073146,0.14852,0.22842,0.31528,0.41176,0.52078,0.64567,0.79024,0.95889,1.1568,1.3899,1.6654,1.9916,2.3786,2.838,3.3839,2.838,2.3786,1.9916,1.6654,1.3899,1.1568,0.95889,0.79024,0.64567,0.52078,0.41176,0.31528,0.22842,0.14852,0.073146};
	VERIFY(CompareDecimals(R, realR, 3));
	VERIFY(CompareDecimals(lags, realLags, 4));

	XCorr(x, R, lags, 'b');
	realR = {0.0045716,0.0092825,0.014276,0.019705,0.025735,0.032549,0.040354,0.04939,0.059931,0.072298,0.086868,0.10409,0.12448,0.14866,0.17737,0.21149,0.17737,0.14866,0.12448,0.10409,0.086868,0.072298,0.059931,0.04939,0.040354,0.032549,0.025735,0.019705,0.014276,0.0092825,0.0045716};
	VERIFY(CompareDecimals(R, realR, 3));
	VERIFY(CompareDecimals(lags, realLags, 4));
	
	XCorr(x, R, lags, 'u');
	realR = {0.073146,0.07426,0.076141,0.078821,0.082351,0.086796,0.092238,0.09878,0.10654,0.11568,0.12635,0.13878,0.1532,0.1699,0.1892,0.21149,0.1892,0.1699,0.1532,0.13878,0.12635,0.11568,0.10654,0.09878,0.092238,0.086796,0.082351,0.078821,0.076141,0.07426,0.073146};
	VERIFY(CompareDecimals(R, realR, 3));
	VERIFY(CompareDecimals(lags, realLags, 4));
	
	XCorr(x, R, lags, 'c');
	realR = {0.073146,0.07426,0.076141,0.078821,0.082351,0.086796,0.092238,0.09878,0.10654,0.11568,0.12635,0.13878,0.1532,0.1699,0.1892,0.21149,0.1892,0.1699,0.1532,0.13878,0.12635,0.11568,0.10654,0.09878,0.092238,0.086796,0.082351,0.078821,0.076141,0.07426,0.073146};
	//VERIFY(CompareDecimals(R, realR, 3));
	UppLog() << "\nBypassed XCorr(x, R, lags, 'c')";
	VERIFY(CompareDecimals(lags, realLags, 4));
	
	{
		UVector<double> xx = {658991.77,660487.490817641,662359.786838351,664705.212961601,667564.567557467,670957.290915256,675043.088016135,679798.951258472,685309.494600102,691699.892981593,698872.987597855,707073.351240186,716098.378235089,726036.273293386,736777.111066518,748184.612559005,760077.315531181,772213.108043024,784169.523954539,795542.954686134,805824.048876876,814157.672892616,819942.886316603,822316.30514584,819863.506030399,812295.169931958,797550.808066033,775467.418470448,745146.908819387,705509.331696135,657745.021352108,601978.811420999,540279.148705618,475777.96240604,411456.251009611,351307.521549351,299304.971698198,258195.623064412,228841.76865112,213672.450221289,209525.927948021,216034.98239422,230497.219435427,250138.993611412,273510.862914356,298043.159261679,322562.446687617,345856.891763955,367258.117349896,386021.405783862,402257.789200772,415742.033688155,426454.827373756,434866.295,440633.45421712,444504.657238416,446509.781242969,446925.19475566,446166.607047645,444303.044963984,441792.209621085,438826.020967055,435684.607377599,432590.41284198,429835.140811589,427547.980613649,425936.265558429,425158.800220282,425189.372238094,426196.916690218,428015.780059881};
		UVector<double> yy = {648785.1275,650255.285,652134.2125,654478.7975,657323.07,660748.1075,664793.1675,669556.5475,675069.305,681432.8125,688649.53,696802.4825,705836.4225,715788.865,726512.3125,737959,749854.125,762018.005,774009.5825,785466.2125,795748.1925,804259.28,810151.4925,812566.3925,810473.855,802847.3425,788557.2025,766722.0375,736395.055,697387.3475,649565.87,594190.5525,532669.54,468172.7475,403705.27,343562.165,290965.52,249204.6625,219717.565,203119.945,198464.19,204059.46,217711.64,236969.75,259762.06,284034.88,308454.07,331801.0075,353215.615,372360.975,388831.5975,402657.31,413772.615,422387.74,428627.8375,432771.5025,435044.7475,435736.8275,435121.315,433475.6775,431090.0925,428220.5025,425150.935,422109.35,419348.41,417053.7425,415412.9225,414544.03,414533.78,415400.2125,417112.9875};

		VectorXd rx, ry;
		Copy(xx, rx);
		Copy(yy, ry);
		
		VectorXd RR, lags;		
		XCorr(rx, ry, RR, lags, 'c');
		realR = {0.012148,0.024273,0.036408,0.048586,0.060842,0.073208,0.085715,0.098391,0.11126,0.12433,0.13762,0.15113,0.16486,0.17879,0.1929,0.20716,0.22152,0.23593,0.25031,0.26458,0.27864,0.29238,0.30567,0.31837,0.33035,0.34148,0.35161,0.36066,0.36856,0.3753,0.38097,0.38569,0.3897,0.39331,0.39687,0.40076,0.40535,0.41096,0.41779,0.42603,0.43571,0.44682,0.45927,0.47293,0.48768,0.50335,0.51983,0.53699,0.55475,0.57304,0.59183,0.61108,0.6308,0.65099,0.67165,0.69278,0.71436,0.73637,0.75874,0.78139,0.80419,0.82701,0.84964,0.87191,0.89358,0.91444,0.93427,0.95287,0.97009,0.98581,0.99997,0.98575,0.96998,0.95272,0.9341,0.91427,0.89345,0.87183,0.84966,0.82715,0.80449,0.78186,0.75942,0.73728,0.71552,0.6942,0.67335,0.65298,0.6331,0.61369,0.59476,0.57631,0.55835,0.54094,0.52412,0.508,0.49268,0.47829,0.46497,0.45285,0.44205,0.43265,0.42465,0.418,0.41254,0.40803,0.40415,0.40053,0.3968,0.39262,0.38766,0.38172,0.37466,0.36641,0.357,0.3465,0.33501,0.32268,0.30963,0.29601,0.28197,0.26762,0.25309,0.23848,0.22385,0.20929,0.19485,0.18057,0.16648,0.15261,0.13896,0.12554,0.11234,0.099351,0.086556,0.07393,0.061446,0.049072,0.036775,0.02452,0.012272};
		VERIFY(CompareDecimals(RR, realR, 4));
	}
}

void TestOthers() {
	{
		UppLog() << "\nVerifying SmoothStep():\n";
		
		for (int order : UVector<int>({3, 5, 7})) {
			VERIFY(-10 == SmoothStep(-11., order, -10., 20., -10., 30.));
			VERIFY(-10 == SmoothStep(-10., order, -10., 20., -10., 30.));
			VERIFY( 10 == SmoothStep(  5., order, -10., 20., -10., 30.));
			VERIFY( 30 == SmoothStep( 20., order, -10., 20., -10., 30.));
			VERIFY( 30 == SmoothStep( 21., order, -10., 20., -10., 30.));
		}
	}
	{
		MultiDimMatrixIndex idx(2, 4);
		idx.RowMajor();
		Matrix<double, 2, 4, RowMajor> d;
		double *p = d.data();
		d(1, 2) = 27;
		VERIFY(p[idx(1, 2)] == 27);
	}
	{
		MultiDimMatrixIndex idx(2, 4);
		MatrixXd d(2, 4);
		double *p = d.data();
		d(1, 2) = 27;
		VERIFY(p[idx(1, 2)] == 27);
	}
	{
		MultiDimMatrix<int> mat(1, 2, 1, 4);	
		int id = 0;
		for (int i3 = 0; i3 < mat.size(3); ++i3)
			for (int i1 = 0; i1 < mat.size(1); ++i1)
				mat(0, i1, 0, i3) = id++;
		id = 0;
		for (int i = 0; i < 8; ++i)
			VERIFY(id++ == *(mat.begin() + i));
		
	}
	{
		MultiDimMatrix<VectorXd> mat(1, 2);
		mat(0, 0) = VectorXd::Constant(20, 2.5);	
		mat(0, 1) = VectorXd::Constant(30, 5.2);	
		VERIFY(!IsNull(mat));
		mat(0, 0)(5) = Null;
		VERIFY(IsNull(mat));
		mat(0, 0) = VectorXd();
		VERIFY(IsNull(mat));
	}
	{
		MultiDimMatrix<MatrixXd> mat(1, 2);
		mat(0, 0) = VectorXd::Constant(20, 2.5);	
		mat(0, 1) = VectorXd::Constant(30, 5.2);	
		VERIFY(IsNum(mat));
		mat(0, 0)(5) = std::numeric_limits<double>::quiet_NaN();
		VERIFY(!IsNum(mat));
		mat(0, 0) = VectorXd();
		VERIFY(!IsNum(mat));
	}
	{
		double x, y, dy, d2y;
		
		x = 3;
		LinearInterpolate(x, 2., 4., 1., 5., y, dy);
		VERIFY(EqualDecimals(y, 3., 10));
		VERIFY(EqualDecimals(dy, 2., 10));
		
		x = 4;
		QuadraticInterpolate(x, 2., 4., 5., 1., 5., 3., y, dy, d2y);
		VERIFY(EqualDecimals(y, 5., 10));
		VERIFY(EqualDecimals(dy,  -0.666666666666, 10));
		VERIFY(EqualDecimals(d2y, -2.666666666666, 10));
	}
	{
		UVector<double> x = {0, 2, 4, 6};
		UVector<double> y = {0, 8, 16, 32};
		
		double a, b;
		LinearRegression(x, y, a, b);
		VERIFY(EqualDecimals(a, 5.2, 10));	
		VERIFY(EqualDecimals(b, -1.6, 10));
	}
	{
		MatrixXd a(4, 4);
		a <<  0,  1,  2,  3,
			 10, 11, 12, 13,
			 20, 21, 22, 23,
			 30, 31, 32, 33;
	
		MatrixXd b = a;
		b.array() += 100.;
		
		{
			MatrixXd A = a;
			
			UppLog() << "\nSource matrix:\n";
			UppLog() << a;
			
			Swap(A, 1, 3);
			
			UppLog() << "\nSwap row and col 1 with 3:\n";
			UppLog() << A << "\n";
			
			VERIFY(A(3, 1) == 13);
			VERIFY(A(1, 3) == 31);
		}
		{
			MatrixXd A = a,
					 B = b;
			
			UppLog() << "\nSource matrices:\n";
			UppLog() << a << "\n";
			UppLog() << b;
			
			Swap(A, B, 1, 3);
			
			UppLog() << "\nSwap row and col 1 of 1st matrix with 3 of 2nd matrix:\n";
			UppLog() << A << "\n";		
			UppLog() << B << "\n";			
			
			VERIFY(A(3, 1) == 13);
			VERIFY(B(1, 3) == 131);
		}
		
	}
}

void TestScatter() {		// Functions found in Scatter
	
	UVector<double> t = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1},
					y = {3.70, 1.00, 4.12, 4.09, 4.76, 3.93, 4.59, 7.37, 6.23, 6.79, 8.29};
	
	UVector<double> yres, yres2;
 	MovingAverage(t, y, 0.4, yres);
	MovingAverage(t, y, 0.4, t, yres2, true);
	
	UVector<double> t2(t.size());
	for (int i = 0; i < t.size(); ++i)
		t2[i] = t[i] + .05;
	
	MovingAverage(t, y, 0.4, t2, yres2, true);
}

void TestHomography() {
	UppLog() << "\nVerifying class Homography. They are transformations that describe the motion between two images, when the camera or the observed object moves.\n";
	
	Pointf p01from(14, 28), p01to(34, 29);
	Pointf p02from(288, 28), p02to(281, 25);
	Pointf p03from(288, 218), p03to(263, 130);
	Pointf p04from(14, 218), p04to(16, 207);
	
	Homography<double> h;
	h.QuadToQuad(p01from, p02from, p03from, p04from, p01to, p02to, p03to, p04to);
	
	UppLog() << "\nIt is tested with itself\n";
					 
	Pointf p01toT = h.Transform(p01from);
	Cout() << "Pfrom: " << p01from << ", Pto: " << p01to << ", Must be: " << p01toT << "\n";
	VERIFY(EqualDecimals(p01to.x, p01toT.x, 6) && EqualDecimals(p01to.y, p01toT.y, 6));

	Pointf p02toT = h.Transform(p02from);
	Cout() << "Pfrom: " << p02from << ", Pto: " << p02to << ", Must be: " << p02toT << "\n";
	VERIFY(EqualDecimals(p02to.x, p02toT.x, 6) && EqualDecimals(p02to.y, p02toT.y, 6));
	
	Pointf p03toT = h.Transform(p03from);
	Cout() << "Pfrom: " << p03from << ", Pto: " << p03to << ", Must be: " << p03toT << "\n";
	VERIFY(EqualDecimals(p03to.x, p03toT.x, 6) && EqualDecimals(p03to.y, p03toT.y, 6));
	
	Pointf p04toT = h.Transform(p04from);
	Cout() << "Pfrom: " << p04from << ", Pto: " << p04to << ", Must be: " << p04toT << "\n";
	VERIFY(EqualDecimals(p04to.x, p04toT.x, 6) && EqualDecimals(p04to.y, p04toT.y, 6));
	
	UppLog() << "\nOther examples\n";
	
	Pointf p1from(85, 150);
	Pointf p1to = h.Transform(p1from);
	Cout() << "Pfrom: " << p1from << ", Pto: " << p1to << "\n";
	
	Pointf p2from(185, 100);
	Pointf p2to = h.Transform(p2from);
	Cout() << "Pfrom: " << p2from << ", Pto: " << p2to << "\n";
	
	Pointf p3from(285, 200);
	Pointf p3to = h.Transform(p3from);
	Cout() << "Pfrom: " << p3from << ", Pto: " << p3to << "\n";
	
	Pointf p4from(185, 250);
	Pointf p4to = h.Transform(p4from);
	Cout() << "Pfrom: " << p4from << ", Pto: " << p4to << "\n";
	
	UppLog() << "\nNow reconstructing an image\n";
	
	Size sz(700, 250);
	 
	ImageBuffer ib(sz);
	for(int y = 0; y < sz.cy; y++)
		Fill(ib[y], White(), sz.cx);
		
	BufferPainter bp(ib);
 	bp.Begin();
  	bp.DrawText(50, 50, "Hello U++!", Arial(120).Bold(), Black());
 	bp.End();
 	
    Image original = ib;
    
    String dir = AppendFileNameX(GetExeFolder(), "STEM4U_Demo");
	RealizeDirectory(dir);
	PNGEncoder().SaveFile(AFX(dir, "Homography_original.png"), original);
	
	Point orig00(0, 0),   orig10(sz.cx-1, 0), orig11(sz.cx-1, sz.cy-1), orig01(0, sz.cy-1);
	Point dest00(0, 200), dest10(600, 0), dest11(700, 500), dest01(20, 400);
	
	Image deformed = ApplyHomography(original, Green(), 
				orig00, orig10, orig11, orig01,
			   	dest00, dest10, dest11, dest01, true);
			   							  
	PNGEncoder().SaveFile(AFX(dir, "Homography_deformed.png"), deformed);	
	
  	Image reconstr = ApplyHomography(deformed, Green(), 
  				dest00, dest10, dest11, dest01,
  				orig00, orig10, orig11, orig01, true);
			   							  
	PNGEncoder().SaveFile(AFX(dir, "Homography_reconstructed.png"), reconstr);
	
	Image reconstr2= ApplyHomography(deformed, Green(), 
  				dest00, dest10, dest11, dest01, sz, true);
			   							  
	PNGEncoder().SaveFile(AFX(dir, "Homography_reconstructed 2.png"), reconstr2);
}

void TestLocalFitting(bool test);
void TestMooring(bool test);
void TestButterworth(bool test);
void TestCombinations();
void TestRational();
void TestPolynomial();
	    	
CONSOLE_APP_MAIN 
{
	StdLogSetup(LOG_COUT|LOG_FILE);
	SetExitCode(0);
	
	UppLog() << "STEM4U demo and test";
	
	try {
		bool test = CommandLine().size() > 0 && CommandLine()[0] == "-test";
		
		TestHomography();
		TestScatter();
		TestVectorMatrixHealing();
		TestRootfinding();
		TestOthers();
		TestCombinations();
		TestXCorr();
		TestMooring(test);
		TestLocalFitting(test);
		TestButterworth(test);
		TestTravellingSalesman();
		TestShortestPath();
		//TestRational();
		TestDAESolver();
		TestIntInf();
	    //TestPolynomial();
	    TestIntegral();
	    TestIntegralSinCos();
	    TestSeaWaves();
	  
	    UppLog() << "\n\nAll tests passed\n";
   	} catch (Exc e) {
		UppLog() << "\nError: " << e << "\n";  
		SetExitCode(1);
	} catch (...) {
		UppLog() << "\nUnknown Error\n";  
		SetExitCode(1);
	}
	 
	#ifdef flagDEBUG
	UppLog() << "\n";
	Cout() << "\nPress enter key to end";
	ReadStdIn();
	#else
	Cout() << "\nPress enter key to end";
	#endif   
}

