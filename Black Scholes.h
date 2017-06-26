//Black Scholes formula for European options

#include <vector>
#include <math.h>
#include <algorithm> 


// approximation of N() 
double N(double x)
{
	const double d1 = 0.0498673470;
	const double d2 = 0.0211410061;
	const double d3 = 0.0032776263;
	const double d4 = 0.0000380036;
	const double d5 = 0.0000488906;
	const double d6 = 0.0000053830;
	double ax = abs(x);
	double n = 1.0 - 1.0 / 2.0*pow((1 + d1*ax + d2*ax*ax + d3*pow(ax, 3) + d4*pow(ax, 4) + d5*pow(ax, 5) + d6*pow(ax, 6)), -16.0);
	if (x<0) n = 1 - n;
	return n;
}

// BS formula for European Call/Put option
double EuropeanBS(double S0, double K, double T,  double r, double sigma, char OpType)
{
	double d = (log(S0 / K) + T*(r + 0.5*sigma*sigma)) / (sigma*sqrt(T));
	double Call = S0*N(d) - exp(-r*T)*K*N(d - sigma*sqrt(T));
	if (OpType == 'C')
		return Call;
	else
		return Call - S0 + K*exp(-r*T);
}
