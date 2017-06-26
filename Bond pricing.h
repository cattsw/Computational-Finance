//Vasicek model for short rate

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <algorithm> 
#include <iostream>
#include "Short rate model.h"


using namespace std; 


// Price of discount bond by monte carlo simulation
double DBond_MC(char model,double F=1000.0, double r0 = 0.05, double x0=0., double y0=0., double rbar = 0.05, double kappa = 0.82, double sigma = 0.1, double T=.5,
			int NPaths=1000, int NSteps = 180)
{
	// model= 'V': Vasicek, model='C': CIR , model='G': G2++

	double dt = .5 / NSteps; // Time increment
	// Initialize the bond price
	double P = 0.0;

	for (int i = 0; i < NPaths; i++)
	{
		vector<double> r;
		switch (model) {
		case 'V': {
			r = RatePath_V(r0, rbar, kappa, sigma, T, int(2 * NSteps*T)); }
			break; 
		case 'C': {
			r = RatePath_CIR(r0, rbar, kappa, sigma, T, int(2 * NSteps*T)); }
			break; 
		case 'G':
			r = RatePath_G2(T, int(NSteps * 2 * T),0.03,r0,x0,y0)[0];
			break;
		default:            // wrong input
			cout << "Error, choose V/C/G \n";
			break;
		}

		double  R = 0.0;// Initialize the integral of rates for each path
		for (int t = 0; t < r.size(); t++)
		{
			R = R + r[t];
		}
		R = R * dt;
		P = P + exp(-R);
	}
	return F*P / NPaths;
}

// Price of coupon paying bond
double CBond_MC( vector<double> Ci = { 30.0,30.0,30.0,30.0,30.0,30.0,30.0,1030.0 },
				double r0 = 0.05, double rbar = 0.05, double kappa = 0.82, double sigma = 0.1, 
	            vector<double> Ti = { 0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0},
				int NPaths = 1000, int NSteps = 180)
{
	// Ci = series of coupon payment
	// Ti= time for each coupon payment
	// NSteps = number of time step for each 0.5 year
	double T = Ti.back();
	double dt = 0.5 / NSteps; // Time increment
	// Initialize the bond price
	double P = 0.0;
	int cn = Ci.size(); //number of coupon payment
	for (int i = 0; i < NPaths; i++)
		{
		double C = 0.0;
		vector<double> r = RatePath_V(r0, rbar, kappa, sigma, T, int(NSteps*2*T));
		for (int t=0; t<cn; t++) 
		{
			double R = 0.0;// Initialize the integral of rates for each path
			for (int j = 0; j <= int(2*Ti[t]*NSteps); j++)
			{
				R = R + r[j];
			}
			R = R * dt;
			C = C + Ci[t]*exp(-R);
		}
		P = P + C;
	}
	return P / NPaths;
}


// Price of European Call on discount bond
double EOpt_DBond_MC(char OptType, char model, double F , double K , double T , double S , 
			double r0=0.05 ,double x0=0.,double y0=0., double rbar=0.05 , double kappa=0.82 , double sigma=0.1 , int NPaths = 1000, int NSteps = 180)
{
	//Initialize option and bond price
	double Opt = 0.;
	double P = 0.;
	double dt = 0.5 / NSteps; // Time increment
	double x=0., y=0.;

	for (int i = 0; i < NPaths; i++)
	{
		// rate path 
		vector<double> r;
		vector< vector<double>> rxy;
		switch (model) {
		case 'V':
			r = RatePath_V(r0, rbar, kappa, sigma, T, int(NSteps * 2 * T));
			break;
		case 'C':
			r = RatePath_CIR(r0, rbar, kappa, sigma, T, int(NSteps * 2 * T));
			break;
		case 'G':
			 rxy = RatePath_G2(T, int(NSteps * 2 * T), 0.03, r0, x0, y0);
			r = rxy[0];
			x = rxy[1].back();
			y = rxy[2].back();
			break;
		default:            // wrong input
			cout << "Error, choose V/C/G \n";
			break;
		}

		double R = 0.0;// Initialize the integral of rates for each path
		for (int j = 0; j <= int(NSteps * 2 * T); j++)
		{
			R = R + r[j];
		}
		R = R * dt;
		// bond price for each path
		P = DBond_MC(model, F, r[int(NSteps * 2 * T)], x,y,rbar, kappa, sigma, S - T);
		//	cout << exp(-R1)*max(P - K, 0.0) << endl;
		if (OptType == 'C')
			Opt = Opt + exp(-R)*max(P - K, 0.0);
		else
			Opt = Opt + exp(-R)*max(K - P, 0.0);
	}
	return Opt / NPaths;
}



// Price of European Call on coupon bond
double ECall_CBond_MC(vector<double> Ci = { 30.0,30.0,30.0,30.0,30.0,30.0,30.0,1030.0 },
					double K = 980.0, double r0 = 0.05, double rbar = 0.05, double kappa = 0.82, double sigma = 0.1,
					double T = 0.25, vector<double> Si = { 0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0 }, int NPaths = 1000, int NSteps = 180)
{
	//Initialize call and bond price
	double Call = 0.;
	double dt = 0.5 / NSteps; // Time increment
	double S = Si.back();
	int cn = Ci.size(); //number of coupon payment
	for (int i = 0; i < NPaths; i++)
	{
		// rate path 
		vector<double> r = RatePath_V(r0, rbar, kappa, sigma, T, int(NSteps * 2 * T));
		double R = 0.0;// Initialize the integral of rates for each path
		for (int j = 0; j <= int(NSteps * 2 * T); j++)
		{
			R = R + r[j];
		}
		R = R * dt;
		// bond price for each path
		vector<double> Ti(cn);
		for (int t = 0; t < cn; t++)
		{
			Ti[t] = Si[t] - T;
		}
		double P = CBond_MC(Ci , r[int(NSteps * 2 * T)],  rbar, kappa , sigma , Ti);
	//	cout << P << endl;
		Call = Call + exp(-R)*max(P - K, 0.0);
	}
	return Call / NPaths;
}

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


// Explicit formula for discounted bond price fom Vasicek model
double B(double t, double T, double kappa = 0.82)
{
	return 1. / kappa*(1. - exp(-kappa*(T - t)));
}

double A(double t, double T, double rbar = 0.05, double kappa = 0.82, double sigma = 0.1)
{
return exp((rbar - sigma*sigma / (2.*kappa*kappa))*(B(t, T, kappa) - (T - t)) - sigma*sigma*B(t, T, kappa)*B(t, T, kappa) / (4.*kappa));
}

double DBond_V(double t, double T, double rt = 0.05, double rbar = 0.05, double kappa = 0.82, double sigma = 0.1)
{
	return A(t, T)*exp(-B(t, T)*rt);
}

//Explicit formula for European Call/Put option on discounted bond
double ECall_DBond_Explicit_V(double t, double T, double S, double K, double rt = 0.05, double rbar = 0.05, double kappa = 0.82, double sigma = 0.1)
{
	double sigmap = sqrt((1. - exp(-2.*kappa*(T - t))) / (2.*kappa)) * (1. - exp(-kappa*(S - T)))/ kappa*sigma;
	double d1 = log(DBond_V(t, S) /( K* DBond_V(t, T))) / sigmap + sigmap / 2.;
	double d2 = d1 - sigmap;
	return  DBond_V(t, S)*N(d1) - K*DBond_V(t, T)*N(d2);
}

// Bisection Algorithm for root search 
double Root(double a, double b, vector<double> Ci = { 30.0,30.0,30.0,30.0,30.0,30.0,30.0,1030.0 }, double K = 980.0,
	 double T = 0.25,vector<double> Si = { 0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0 })
{
	const int n = 10000000;  //number of max iterations
	double epsilon = 0.001;   // Tolerance
	double mid = 0.5 * (a + b);
	int cn = Ci.size();
	double midK = 0.,lK=0.,hK=0.;
	for (int i = 0; i < cn; i++) {
	    midK= midK+Ci[i]* DBond_V(T, Si[i],mid);
		lK = lK+ Ci[i] * DBond_V(T, Si[i], a);
		hK = hK+Ci[i] * DBond_V(T, Si[i], b);
	}
	// make sure choose the appropriate interval
	double  l = K -lK;
	double h =K -hK;
	if (l*h > 0)
	{
		printf("wrong interval\n");
		return -1;
	}

	// While the difference between midP and the market price is greater than epsilon, keep subdividing the interval 
	int i = 0;
	do {
		if (midK < K) {
			a = mid;
		}
		else
		{
			b = mid;
		}
		mid = 0.5 * (a + b);
		midK = 0.;
		for (int i = 0; i < cn; i++) {
			midK = midK+Ci[i] * DBond_V(T, Si[i], mid);
		}
		i = i + 1;

	} while ((abs(midK - K) > epsilon) & (i <= n));

	return mid;
}

double ECall_CBond_Explicit_V(vector<double> Ci = { 30.0,30.0,30.0,30.0,30.0,30.0,30.0,1030.0 }, double K = 980.0,
	double r0 = 0.05, double rbar = 0.05, double kappa = 0.82, double sigma = 0.1,double T=0.25,
	vector<double> Si = { 0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0 })
{
	int cn = Ci.size();
	double Call =0.;
	double r_star = Root(1., 0.);
	//cout << r_star << endl;
	for (int i = 0; i < cn; i++)
	{
		double Ki = DBond_V(T, Si[i], r_star);
	//	cout << ECall_DBond_V(0., T, Si[i], Ki) << endl;
		Call = Call + Ci[i]*ECall_DBond_Explicit_V(0., T, Si[i], Ki);
	}
	return Call;
}



// Explicit formula for discounted bond price fom CIR model
double DBond_CIR(double t, double T, double rt = 0.05, double rbar = 0.055, double kappa = 0.92, double sigma = 0.12)
// this is the original CIR formulation of their term structure.
{
	double sigma_sqr = pow(sigma, 2);
	double gamma = sqrt(pow(kappa, 2) + 2.0*sigma_sqr);
	double denum = (gamma + kappa)*(exp(gamma*(T - t)) - 1) + 2 * gamma;
	double p = 2 * kappa*rbar / sigma_sqr;
	double enum1 = 2 * gamma*exp(0.5*(kappa + gamma)*(T - t));
	double A = pow((enum1 / denum), p);
	double B = (2 * (exp(gamma*(T - t)) - 1)) / denum;
	double DBond = A*exp(-B*rt);
	return DBond;
};


//Implicit Finite difference for bond option

double ECall_EFD(double K, double T = 0.5, double r0 = 0.05, double rbar = 0.055, double kappa = 0.92, double sigma = 0.12, char CorP = 'C', int n = 180, int m = 180)
{
	int NPaths = m;				// Number of stock price paths
	double dt = T / n;			// Time increment
	int NSteps = n;		  // Number of time steps
						  //Transformed PDE 
	double mu;	  // Drift for r process
	double dr = sigma * sqrt(dt) ;

	double pu;      //Up probability
	double pm;      //Middle probability
	double pd;       //Down probability

	 //Initialize r and option values.
	vector<double> r(2 * NPaths + 1, 0.0);
	vector<vector<double>> Opt(2 * NPaths + 1, vector<double>(NSteps + 1, 0.0));

	// Indices for r step
	vector<double> J(2 * NPaths + 1);
	for (int i = NPaths; i >= -NPaths; --i)
		J[NPaths - i] = i;

	// r at maturity
	for (int i = 0; i < r.size(); ++i)
		r[i] = r0 * exp(J[i] * dr);

	// Option price at maturity
	for (int i = 0; i < Opt.size(); i++)
		Opt[i][Opt[0].size() - 1] = ((CorP == 'P') ? max(K -r[i], 0.0) : max( r[i] - K, 0.0));

	//Initialize required matrices for solving the triangular system of equations.
	vector <vector<double>> pmp(2 * NPaths, vector<double>(NSteps));
	vector <vector<double>> pp(2 * NPaths, vector<double>(NSteps));
	vector<double> C(2 * NPaths + 1);

	//Upper and lower bounds 
	double lambda_L = (CorP == 'P') ? max(K -  r[2*NPaths-1], 0.0) : 0.0;
	double lambda_U = (CorP == 'C') ? max( r[0] - K, 0.0) : 0.0;

	//Solve the triangular system of equations.
	for (int j = NSteps - 1; j >= 0; --j)
	{
		for (int i = 0; i < 2 * NPaths + 1; ++i)
		{
			if (i == 0)
				C[i] = lambda_U;
			else if (i < 2 * NPaths)
				C[i] = Opt[i][j + 1];
			else
				C[i] = lambda_L;
		}
		mu = kappa*(rbar - r[2 * NPaths - 1]) ;	  // Drift for r process
		pm = 1.0 + dt * (sigma * sigma*r[2 * NPaths - 1]) / (dr * dr) + r[2 * NPaths - 1] * dt;      //Middle probability
		pd = -0.5 * dt* ((sigma * sigma*r[2 * NPaths - 1]) / (dr * dr) - mu / dr);       //Down probability

		pmp[2 * NPaths - 1][j] = pd + pm;
		pp[2 * NPaths - 1][j] = C[2 * NPaths - 1] + pd*lambda_L;
		for (int i = 2 * NPaths - 2; i > 0; --i)
		{
			mu = kappa*(rbar - r[i]) ;	  // Drift for stock process
			pu = -0.5 * dt * ((sigma * sigma*r[i]) / (dr * dr) + mu / dr);      //Up probability
			pm = 1.0 + dt * (sigma * sigma*r[i]) / (dr * dr) + r[i] * dt;      //Middle probability
			pd = -0.5 * dt* ((sigma * sigma*r[i]) / (dr * dr) - mu / dr);       //Down probability

			pmp[i][j] = pm - pu*pd / (pmp[i + 1][j]);
			pp[i][j] = C[i] - pp[i + 1][j] * pd / pmp[i + 1][j];
		}

		mu = kappa*(rbar - r[0]) ;
		pu = -0.5 * dt * ((sigma * sigma*r[0]) / (dr * dr) + mu / dr);
		Opt[0][j] = (pp[1][j] + pmp[1][j] * lambda_U) / (pu + pmp[1][j]);
		Opt[1][j] = Opt[0][j] - lambda_U;

		for (int i = 2; i < 2 * NPaths; i++)
		{
			mu = kappa*(rbar - r[i]) ;
			pu = -0.5 * dt * ((sigma * sigma*r[i]) / (dr * dr) + mu / dr);
			Opt[i][j] = (pp[i][j] - pu*Opt[i - 1][j]) / pmp[i][j];
		}

		Opt[2 * NPaths][j] = Opt[2 * NPaths - 1][j] - lambda_L;


	}
	return Opt[NPaths][0];
}




// Explicit formula for discounted bond price fom g2++ model

double DBond_G2(double t, double T, double phi = 0.03, double r0 = 0.03, double a = 0.1, double b = 0.3,
	double sigma1 = 0.03, double sigma2 = 0.08, double pho = 0.7)
{
	double V = sigma1*sigma1 / (a*a)*(T - t + 2.*exp(-a*(T - t)) / a - exp(-2.*a*(T - t)) / (2.* a) - 1.5 / a) +
		sigma2*sigma2 / (b*b)*(T - t + 2.*exp(-b*(T - t)) / b - exp(-2.*b*(T - t)) / (2.* b) - 1.5 / b) +
		2 * pho*sigma1*sigma2 / (a*b)*(T - t + (exp(-a*(T - t)) - 1) / a + (exp(-b*(T - t)) - 1) / b - (exp(-(a + b)*(T - t)) - 1) / (a + b));

	return exp(-phi*(T - t) + 0.5*V);
}

//Explicit formula for European Call/Put option on discounted bond by G2 model
double EPut_DBond_Explicit_G2(double t, double T, double S, double L, double K, double phi = 0.03, double r0 = 0.03, double a = 0.1, double b = 0.3,
				double sigma1 = 0.03, double sigma2 = 0.08, double pho = 0.7)
{
	double Sigma =sqrt( sigma1*sigma1 / (2 * a*a*a)*(1 - exp(-a*(S - T)))*(1 - exp(-a*(S - T)))*(1 - exp(-2.*a*(T - t))) +
		sigma2*sigma2 / (2 * b*b*b)*(1 - exp(-b*(S - T)))*(1 - exp(-b*(S - T)))*(1 - exp(-2.*b*(T - t))) +
		2 * pho*sigma1*sigma2 / (a*b*(a + b))*(1 - exp(-a*(S - T)))*(1 - exp(-b*(S - T)))*(1 - exp(-(a + b)*(T - t))));
	double d1 = log(L*DBond_G2(t, S) / (K*DBond_G2(t, T))) / Sigma + Sigma / 2.;
	double d2 = d1 - Sigma;
	return -L*DBond_G2(t, S)*N(-d1) +DBond_G2(t, T)*K*N(-d2);
}



