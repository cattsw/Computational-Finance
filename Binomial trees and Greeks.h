

#include <vector>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <random>
#include <algorithm> 
using namespace std;


// Q1:Binomial trees with different parameters

double BinomialCall(double S0, double K, double T, double r, double sigma, int n, char ParType)
{
	// Tree Parameters	
	double dt = T / n;
	double u, d, p;

	// Switch parameters for difference cases
	switch (ParType) {
	case 'a':
	{
		double c = 1.0 / 2 * (exp(-r*dt) + exp((r + sigma*sigma)*dt));
		d = c - sqrt(c*c - 1);
		u = 1 / d;
		p = (exp(r*dt) - d) / (u - d);
	}
		break;
	case 'b':
	{
		u = exp(r*dt)*(1.0 + sqrt(exp(sigma*sigma*dt) - 1.0));
		d = exp(r*dt)*(1.0 - sqrt(exp(sigma*sigma*dt) - 1.0));
		p = 1.0 / 2;
	}
		break;
	case 'c': //JR model
	{
		u = exp((r - 1.0 / 2 * sigma*sigma)*dt + sigma*sqrt(dt));
		d = exp((r - 1.0 / 2 * sigma*sigma)*dt - sigma*sqrt(dt));
		p = 1.0 / 2;
	}
		break;
	case 'd': // CRR model
	{
		u = exp(sigma*sqrt(dt));
		d = exp(-sigma*sqrt(dt));
		p = 1.0 / 2 + 1.0 / 2 * ((r - 1.0 / 2 * sigma*sigma)*sqrt(dt) / sigma);
	}
		break;
	}

	// Stock price paths
	vector<vector<double> >  S(n + 1, vector<double>(n + 1));
	S[0][0] = S0;
	for (int i = 0; i <= n; i++)
	{
		for (int j = i; j <= n; j++)
		{
			S[i][j] = S0*pow(u, j - i)*pow(d, i);
		}
	}
	//  Price paths for European Call
	vector<vector<double>> CallOp(n + 1, vector<double>(n + 1));
	for (int i = 0; i <= n; i++)
	{
		CallOp[i][n] = max(S[i][n] - K, 0.);
	}
	for (int j = n - 1; j >= 0; j--) 
	{
		for (int i = 0; i <= j; i++)
		{
			CallOp[i][j] = exp(-r*dt)*(p*CallOp[i ][j + 1] + (1.0 - p)*CallOp[i + 1][j+1]);
		}
	}

	return CallOp[0][0];
}


// Q2:Implied volatility by Bisection Algorithm

//use the approximation of N() described in chapter 3
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

// BS formula for European Call option
double EuropeanCallBS(double S_0,  double X, double r, double T, double Sigma)
{
	double d = (log(S_0 / X) + T*(r + 0.5*Sigma*Sigma)) / (Sigma*sqrt(T));
	double C = S_0*N(d) - exp(-r*T)*X*N(d - Sigma*sqrt(T));
	return  C;
}

// Bisection Algorithm for root search (implied volatility from market price)
double BisecIV(double S, double K, double r, double T,
	double a, double b, double MktPrice) 
{
	const int n = 1000000;  //number of max iterations
	double epsilon=0.005;   // Tolerance
	double mid = 0.5 * (a+b);
	double midP = EuropeanCallBS(S, K, r, T, mid);

	// make sure choose the appropriate interval
	double  l = MktPrice - EuropeanCallBS(S, K, r, T, a); 
	double h= MktPrice - EuropeanCallBS(S, K, r, T, b); 
	if (l*h > 0)
		return -1;

	// While the difference between midP and the market price is greater than epsilon, keep subdividing the interval 
	int i = 0;
	do {
		if (midP < MktPrice) {
			a = mid;
		}
		else
		{
			b = mid;
		}
		mid = 0.5 * (a + b);
		midP = EuropeanCallBS(S, K, r, T, mid);
		i = i + 1;

	} while ((abs(midP - MktPrice) > epsilon) & (i<=n));

	return mid;
}


//Q3: binomial method (choose CRR model) to estimate Greeks
vector<double> BinomialGreeks(double S0, double K, double T, double r, double sigma, int n)
{
	// Tree Parameters	
	double dt = T / n;
	double u = exp(sigma*sqrt(dt));
	double d = exp(-sigma*sqrt(dt));
	double p = 1.0 / 2 + 1.0 / 2 * ((r - 1.0 / 2 * sigma*sigma)*sqrt(dt) / sigma);
	vector<double> greeks(5, 0.0); 
	
	// Stock price paths
	vector<vector<double> >  S(n + 1, vector<double>(n + 1));
	S[0][0] = S0;
	for (int i = 0; i <= n; i++)
	{
		for (int j = i; j <= n; j++)
		{
			S[i][j] = S0*pow(u, j - i)*pow(d, i);
		}
	}
	//  Price paths for European Call
	vector<vector<double>> CallOp(n + 1, vector<double>(n + 1));
	for (int i = 0; i <= n; i++)
	{
		CallOp[i][n] = max(S[i][n] - K, 0.);
	}
	for (int j = n - 1; j >= 0; j--)
	{
		for (int i = 0; i <= j; i++)
		{
			CallOp[i][j] = exp(-r*dt)*(p*CallOp[i][j + 1] + (1.0 - p)*CallOp[i + 1][j + 1]);
		}
	}

	//Calculate greeks 
	if (T > 0)
	{
		double delta = (CallOp[0][1] - CallOp[1][1]) / (S[0][1] - S[1][1]);
		double deltaUp = (CallOp[0][2] - CallOp[1][2]) / (S[0][2] - S[1][2]);
		double deltaDown = (CallOp[1][2] - CallOp[2][2]) / (S[1][2] - S[2][2]);
		double gamma = (deltaUp - deltaDown) / ((S[0][2] - S[2][2])/2.0);
		double theta = (CallOp[0][0] - CallOp[1][2]) / (2 * dt);

		// Calculate vega and rho with re -evaluation
		double uv = exp((sigma + 0.001)*sqrt(dt));
		double dv = exp(-(sigma + 0.001)*sqrt(dt));
		double pv = 1.0 / 2 + 1.0 / 2 * ((r - 1.0 / 2 * (sigma + 0.001)*(sigma + 0.001))*sqrt(dt) / (sigma + 0.001));
		double pr = 1.0 / 2 + 1.0 / 2 * ((r + 0.001 - 1.0 / 2 * sigma*sigma)*sqrt(dt) / sigma);
		vector<vector<double> >  Sv(n + 1, vector<double>(n + 1));
		Sv[0][0] = S0;
		for (int i = 0; i <= n; i++)
		{
			for (int j = i; j <= n; j++)
			{
				Sv[i][j] = S0*pow(uv, j - i)*pow(dv, i);
			}
		}
		vector<vector<double>> CallOpv(n + 1, vector<double>(n + 1));
		vector<vector<double>> CallOpr(n + 1, vector<double>(n + 1));
		for (int i = 0; i <= n; i++)
		{
			CallOpv[i][n] = max(Sv[i][n] - K, 0.);
			CallOpr[i][n] = max(S[i][n] - K, 0.);
		}
		for (int j = n - 1; j >= 0; j--)
		{
			for (int i = 0; i <= j; i++)
			{
				CallOpv[i][j] = exp(-r*dt)*(pv*CallOpv[i][j + 1] + (1.0 - pv)*CallOpv[i + 1][j + 1]);
				CallOpr[i][j] = exp(-(r + 0.001)*dt)*(pr*CallOpr[i][j + 1] + (1.0 - pr)*CallOpr[i + 1][j + 1]);
			}
		}

		double vega = (CallOpv[0][0] - CallOp[0][0]) / 0.001;
		double rho = (CallOpr[0][0] - CallOp[0][0]) / 0.001;
		greeks[0] = delta;
		greeks[1] = gamma;
		greeks[2] = theta / 365.0;
		greeks[3] = vega / 100.0;
		greeks[4] = rho / 100.0;
	}
	return greeks;
}


//Q4: Binomial Tree for both European and American options

double Binomial(double S0, double K, double T, double r, double sigma, int n, char CorP, char EorA)
{
	// Tree Parameters	
	double dt = T / n;
	double u = exp(sigma*sqrt(dt));
	double d = exp(-sigma*sqrt(dt));
	double p = (exp(r*dt) - d) / (u - d);
	
	// Stock price paths
	vector<vector<double> >  S(n + 1, vector<double>(n + 1));
	S[0][0] = S0;
	for (int i = 0; i <= n; i++)
	{
		for (int j = i; j <= n; j++)
		{
			S[i][j] = S0*pow(u, j - i)*pow(d, i);
		}
	}
	//  Price paths for option with different styles
	vector<vector<double>> Opt(n + 1, vector<double>(n + 1));
	for (int i = 0; i <= n; i++)
	{
		switch (CorP) {
		case 'C': Opt[i][n] = max(S[i][n] - K, 0.0); break;
		case 'P': Opt[i][n] = max(K - S[i][n], 0.0); break;
		}
	}
	for (int j = n - 1; j >= 0; j--)
	{
		for (int i = 0; i <= j; i++)
		{
			switch (EorA) {
			case 'A':
				if (EorA == 'C') Opt[i][j] = max(S[i][j] - K, exp(-r*dt)*(p*Opt[i][j + 1] + (1.0 - p)*Opt[i + 1][j + 1]));
				else Opt[i][j] = max(K - S[i][j], exp(-r*dt)*(p*Opt[i][j + 1] + (1.0 - p)*Opt[i + 1][j + 1]));
				break;
			case 'E': Opt[i][j] = exp(-r*dt)*(p*Opt[i][j + 1] + (1.0 - p)*Opt[i + 1][j + 1]); break;
			}
			
		}
	}

	return Opt[0][0];
}

//Q5 Trinomial tree for European Call option
double TrinomialCall(double S0, double K, double T, double r, double sigma, int n, char ParType)
{
	// Tree Parameters	
	double dt = T / n;
	double pd, pu, pm; 
	// Stock price paths
	vector<vector<double> >  S(2 * n + 1, vector<double>(n + 1));
	//  Price paths for European Call
	vector<vector<double>> CallOp(2 * n + 1, vector<double>(n + 1));

	// Switch parameters for difference cases
	switch (ParType) {
	case 'a':
	{
		double d = exp(-sigma*sqrt(3 * dt));
		double u = 1 / d;
		pd = (r*dt*(1 - u) + r*r*dt*dt + sigma*sigma*dt) / (u - d) / (1 - d);
		pu = (r*dt*(1 - d) + r*r*dt*dt + sigma*sigma*dt) / (u - d) / (u - 1);
		pm = 1 - pu - pd;
		// Stock price paths
		for (int j = 0; j <= n; j++)
		{
			for (int i = 0; i <= 2 * j; i++)
			{
				S[i][j] = S0*pow(u, j - i);
			}
		}
		for (int i = 0; i <= 2 * n; i++)
		{
			CallOp[i][n] = max(S[i][n] - K, 0.0);
		}
	}
	break;
	case 'b':
	{
		double dx = sigma*sqrt(3 * dt);
		pd = 1.0 / 2 * ((sigma*sigma*dt + (r - sigma*sigma / 2.0)*(r - sigma*sigma / 2.0)*dt*dt) / (dx*dx) - (r - sigma*sigma / 2.0)*dt / dx);
		pu = 1.0 / 2 * ((sigma*sigma*dt + (r - sigma*sigma / 2.0)*(r - sigma*sigma / 2.0)*dt*dt) / (dx*dx) + (r - sigma*sigma / 2.0)*dt / dx);
		pm = 1 - pu - pd;
		// Stock price paths
		for (int j = 0; j <= n; j++)
		{
			for (int i = 0; i <= 2 * j; i++)
			{
				S[i][j] = log(S0) + (j - i)*dx;
			}
		}
		for (int i = 0; i <= 2 * n; i++)
		{
			CallOp[i][n] = max(exp(S[i][n] )- K, 0.0);
		}
	}
	break;
	}
	
	// Backward price paths
	
	for (int j = n - 1; j >= 0; j--)
	{
		for (int i = 0; i <= 2*j; i++)
		{
			CallOp[i][j] = exp(-r*dt)*(pu*CallOp[i][j + 1] + pm*CallOp[i+1][j + 1]+pd*CallOp[i + 2][j + 1]);
		}
	}

	return CallOp[0][0];
	}

//Q6: Halton's low discrepancy sequences to price Euopean call option
// Halton sequences: input:  number of points, base
vector<double> HaltonSeq(int num, int base)
{
	vector<double> hs;
	double b_inv = 1.0 / base;
	int k = 1;
	//  find the set of coefficients 
	int r_max = floor(log(k) / log(base));
	vector<int> a(r_max + 1, 0);
	int am = pow(base, r_max);
	for (int j = 0; j <= r_max; j++) {
		a[j] = floor(k / am);
		k = k - am * a[j];
		am = am / base;
	}

	for (int i = 0; i < num; i++)
	{
		int nr = a.size();
		b_inv = 1.0 / base;
		double hn = 0.0;
		for (int j = 0; j < nr; j++) {
			hn += b_inv * a[j];
			b_inv /= base;
		}
		hs.push_back(hn);

		//increment the base-b expansion by one
		bool ok = true;
		for (int i = 0; i < nr; i++) {
			if (ok) {
				if (a[i] == base - 1)
					a[i] = 0;
				else {
					a[i] += 1;
					ok = false;
				}
			}
		}
		if (ok) a.push_back(1);
	}
	return hs;
}

// Generate normally distributed random numbers using Halton's sequence
vector<double> Random_Normal_Halton(int N,int b1,int b2)
{
	const double pi = 3.14159265358979323846;
	vector<double> u1(N), u2(N), z(2 * N);
		u1 = HaltonSeq(N, b1);
		u2 = HaltonSeq(N, b2);
		for (int i = 0; i < N; i++)
		{
			z[2 * i] = sqrt(-2.0 * log(u1[i])) * cos(2.0 * pi * u2[i]);
			z[2 * i + 1] = sqrt(-2.0 * log(u1[i])) * sin(2.0 * pi * u2[i]);
		}
	return  z;
}


//price of European Call by Monte Carlo simulation using Halton's sequence
double EuropeanCall_Halton( double S_0,  double K, double T, double r, double Sigma,int N,int b1,int b2)
{
	
	vector<double> ST(2*N, 0.0);		 // Initialize terminal prices S(T)
	vector<double> CallPayoff(2 * N, 0.0);	// Initialize call payoff
	vector<double> z = Random_Normal_Halton(N,  b1,  b2);
	double C = 0.0; // call price
	for (int i = 0; i < 2*N; i++) {
		ST[i] = S_0*exp((r - 0.5*Sigma*Sigma)*T + Sigma*sqrt(T)*z[i]);
		CallPayoff[i] = max(ST[i] - K, 0.0);
		// Simulated prices as discounted average of terminal prices
		C = C + exp(-r*T)*CallPayoff[i];
	}
	
	return C/(2*N);
}

//end