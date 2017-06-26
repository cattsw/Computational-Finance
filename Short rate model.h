//CIR model for short rate

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <algorithm> 

using namespace std;

// Simulate rate paths by Vasicek model
vector<double> RatePath_V(double r0 = 0.05, double rbar = 0.05, double kappa = 0.82, double sigma = 0.1, double T = .5, int NSteps = 180)
{
	double dt = T / NSteps; // Time increment
	vector<double> r(NSteps + 1, 0.0);
	mt19937_64 rng;   // Random number generator 		  
	rng.seed(random_device{}());
	normal_distribution<double> N(0.0, 1.0);  	// Define the standard normal distributions

	r[0] = r0;
	for (int i = 1; i <= NSteps; i++)
	{
		// Euler discretization of r
		double z = sigma*sqrt(dt)*N(rng);
		r[i] = r[i - 1] + kappa*(rbar - r[i - 1])*dt + z;
	}
	return r;
}

// Simulate rate paths by CIR model
vector<double> RatePath_CIR(double r0 = 0.05, double rbar = 0.055, double kappa = 0.92, double sigma = 0.12, double T = 1., int NSteps = 360)
{
	double dt = T / NSteps; // Time increment
	vector<double> r(NSteps + 1, 0.0);
	mt19937_64 rng;   // Random number generator 		  
	rng.seed(random_device{}());
	normal_distribution<double> N(0.0, 1.0);  	// Define the standard normal distributions

	r[0] = r0;
	for (int i = 1; i <= NSteps; i++)
	{
		// Euler discretization of r
		double z = sigma*sqrt(dt)*N(rng);
		r[i] = r[i - 1] + kappa*(rbar - r[i - 1])*dt + sqrt(r[i - 1])*z;
	}
	return r;
}

// Simulate rate paths by G2++ model
vector <vector<double>> RatePath_G2(double T = 1., int NSteps = 360, double phi=0.03, double r0=0.03, double x0=0.,double y0=0., 
			double a= 0.1,double b=0.3,double sigma1 = 0.03, double sigma2=0.08,double pho=0.7)
{
	double dt = T / NSteps; // Time increment
	vector<double> r(NSteps + 1, 0.0);
	vector<double> x(NSteps + 1, 0.0);
	vector<double> y(NSteps + 1, 0.0);
	mt19937_64 rng;   // Random number generator 		  
	rng.seed(random_device{}());
	normal_distribution<double> N1(0.0, 1.0);  	// Define the independent standard normal distributions
	normal_distribution<double> N2(0.0, 1.0);  	

	r[0] = r0;
	x[0] = x0;
	y[0] = y0;
	for (int i = 1; i <= NSteps; i++)
	{
		// Euler discretization of r
		double z1 = N1(rng);
		double z2= N2(rng);
		double x2= pho*z1 + sqrt(1.0 - pho*pho)*z2;
		x[i] = x[i - 1] - a*x[i - 1] * dt + sigma1*sqrt(dt)*z1;
		y[i] = y[i - 1] - b*y[i - 1] * dt + sigma2*sqrt(dt)*x2;
		r[i] = x[i]+y[i]+phi;
	}
	vector <vector<double>> rxy(3, vector<double>(NSteps + 1, 0.0));
	rxy[0] = r; rxy[1] = x; rxy[2] = y;
	return rxy;
}

