//Finite Difference using (transformed) BS PDE 


#include <vector>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <algorithm> 

using namespace std;

//Explicit Finite Difference
double Explicit_FD_X(double S0, double K, double T, double r, double sigma, char CorP, char AorE,int n=250,int m=100)
{
	// CorP: 'C' for Call option, 'P' for Put option
	// AorE: 'A' for American option, 'E' for European option

	int NPaths = m;				// Number of stock price paths
	int NSteps = n;		  // Number of time steps
	double dt = 0.002;
 //Transformed PDE 
	double mu = r - (sigma*sigma) / 2.0;	  // Drift for stock process
	double dx = 100.0*sigma * sqrt(dt) / m;		// Increment for stock price

	double pu = dt*(sigma*sigma / (2.0*dx*dx) + mu / (2.0* dx));         //Up probability
	double pm = 1 - dt*sigma*sigma / (dx*dx) - r*dt;                  //Middle probability
	double pd = dt*(sigma*sigma / (2.0*dx*dx) - mu / (2.0 * dx));     //Down probability

	 //Initialize stock price and option values.
	vector<double> S(2 * NPaths + 1, 0.0);
	vector<vector<double>> Opt(2 * NPaths + 1, vector<double>(NSteps + 1, 0.0));

	// Indices for stock price step
	vector<double> J(2 * NPaths + 1);
	for (int i = NPaths; i >= -NPaths; --i)
		J[NPaths - i] = i;

	// Stock price at maturity
	for (int i = 0; i < S.size(); ++i)
		S[i] = S0 * exp(J[i] * dx);
	
	// Option price at maturity
	for (int i = 0; i < Opt.size(); i++)
		Opt[i][Opt[0].size() - 1] = ((CorP == 'P') ? max(K - S[i], 0.0) : max(S[i] - K, 0.0));

	// Work backwards and obtain matrix of Opt
	for (int j = NSteps-1; j >=0; j--)
	{
		for (int i = 1; i < 2 * NPaths; ++i)
		{
			Opt[i][j] = pu*Opt[i - 1][j + 1] + pm*Opt[i][j + 1] + pd*Opt[i + 1][j + 1];
		}

		if (CorP == 'C') {
			//Lower boundary
			Opt[2 * NPaths ][j] = Opt[2 * NPaths-1][j];
			//Upper boundary
			Opt[0][j] = Opt[1][j] + (S[1] - S[0]);
		}
		else if (CorP == 'P')
		{
			//Lower boundary
			Opt[2 * NPaths ][j] = Opt[2 * NPaths-1][j] + (S[2 * NPaths-1] - S[2 * NPaths]);
			//Upper boundary
			Opt[0][j] = Opt[1][j];
		}
		
		if (AorE=='A'){
		for (int i = 0;i<=2 * NPaths ;++i)
			Opt[i][ j] = ((CorP == 'P') ? max(K - S[i], Opt[i][ j]) : max(S[i]-K, Opt[i][j]));
	     }
		 
	}
	return Opt[NPaths][0];	
}


// Implicit Finite Difference
double  Implicit_FD_X(double S0, double K, double T, double r, double sigma, char CorP, char AorE, int n = 250, int m = 100)
{
	int NPaths = m;				// Number of stock price paths
	double dt = 0.002;			// Time increment
	int NSteps = n;		  // Number of time steps
									  //Transformed PDE 
	double mu = r - (sigma*sigma) / 2;	  // Drift for stock process
	double dx = 100.0*sigma * sqrt( dt)/m;		// Increment for stock price

	double pu = -0.5 * dt * ((sigma * sigma) / (dx * dx) + mu / dx);      //Up probability
	double pm = 1.0 + dt * (sigma * sigma) / (dx * dx) + r*dt;      //Middle probability
	double pd = -0.5 * dt*((sigma * sigma) / (dx * dx) - mu / dx);      //Down probability

																		 //Initialize stock price and option values.
	vector<double> S(2 * NPaths + 1, 0.0);
	vector<vector<double>> Opt(2 * NPaths + 1, vector<double>(NSteps + 1, 0.0));

	// Indices for stock price step
	vector<double> J(2 * NPaths + 1);
	for (int i = NPaths; i >= -NPaths; --i)
		J[NPaths - i] = i;

	// Stock price at maturity
	for (int i = 0; i < S.size(); ++i)
		S[i] = S0 * exp(J[i] * dx);

	// Option price at maturity
	for (int i = 0; i < Opt.size(); i++)
		Opt[i][Opt[0].size() - 1] = ((CorP == 'P') ? max(K - S[i], 0.0) : max(S[i] - K, 0.0));

	//Initialize required matrices for solving the triangular system of equations.
	vector <vector<double>> pmp(2 * NPaths , vector<double>(NSteps ));
	vector <vector<double>> pp(2 * NPaths , vector<double>(NSteps ));
	vector<double> C(2 * NPaths + 1);

	//Upper and lower bounds 
	double lambda_L = (CorP == 'P') ? max(K - S[2 * NPaths], 0.0) : 0.0;
	double lambda_U = (CorP == 'C') ? max(S[0] - K, 0.0) : 0.0;

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

		pmp[2 * NPaths - 1][j] = pd + pm;
		pp[2 * NPaths - 1][j] = C[2 * NPaths - 1] + pd*lambda_L;
		for (int i = 2 * NPaths - 2; i > 0; --i)
		{
			pmp[i][j] = pm - pu*pd / (pmp[i + 1][j]);
			pp[i][j] = C[i] - pp[i + 1][j] * pd / pmp[i + 1][j];
		}

		Opt[0][j] = (pp[1][j] + pmp[1][j] * lambda_U) / (pu + pmp[1][j]);
		Opt[1][j] = Opt[0][j] - lambda_U;

		for (int i = 2; i < 2 * NPaths ; i++)
			Opt[i][j] = (pp[i][j] - pu*Opt[i - 1][j]) / pmp[i][j];

		Opt[2 * NPaths][j] = Opt[2 * NPaths - 1][j] - lambda_L;

		if (AorE == 'A') {
			for (int i = 0; i <= 2 * NPaths; ++i)
				Opt[i][j] = ((CorP == 'P') ? max(K - S[i], Opt[i][j]) : max(S[i] - K, Opt[i][j]));
		}

	}
	return Opt[NPaths][0];
}



//Crank Nicolson
double CrankNicolson_FD_X(double S0, double K, double T, double r, double sigma, char CorP, char AorE, int n = 250, int m = 100)

{
	int NPaths = m;				// Number of stock price paths
	double dt = 0.002;			// Time increment
	int NSteps = n;		  // Number of time steps
	 //Transformed PDE 
	double mu = r  - (sigma*sigma) / 2;	  // Drift for stock process
	double dx = 100.0*sigma * sqrt(dt) / m;			// Increment for stock price
					
	double pu = -0.25 * dt * ((sigma * sigma) / (dx * dx) + mu / dx);      //Up probability
	double pm = 1.0 + dt * (sigma * sigma) / 2.0 / (dx * dx) + r*dt / 2.0;      //Middle probability
	double pd = -0.25 * dt*((sigma * sigma) / (dx * dx) - mu / dx);      //Down probability

	//Initialize stock price and option values.
	vector<double> S(2 * NPaths + 1, 0.0);
	vector<vector<double>> Opt(2 * NPaths + 1, vector<double>(NSteps + 1,0.0));

	// Indices for stock price step
	vector<double> J(2 * NPaths + 1);
	for (int i = NPaths; i >= -NPaths; --i)
		J[NPaths - i] = i;

	// Stock price at maturity
	for (int i = 0; i < S.size(); ++i)
		S[i] = S0 * exp(J[i] * dx);

	// Option price at maturity
	for (int i = 0; i < Opt.size(); i++)
		Opt[i][Opt[0].size() - 1] = ((CorP == 'P') ? max(K - S[i], 0.0) : max(S[i] - K, 0.0)); 

   //Initialize required matrices for solving the triangular system of equations.
	vector <vector<double>> pmp(2 * NPaths + 1, vector<double>(NSteps + 1));
	vector <vector<double>> pp(2 * NPaths + 1, vector<double>(NSteps + 1));
	vector <vector<double>> C(2 * NPaths + 1, vector<double>(NSteps + 1));

	//Solve the triangular system of equations.
	for (int j = NSteps; j > 0; --j)
	{
		pmp[2 * NPaths - 1][j] = pd + pm;
		for (int i = 2 * NPaths - 2; i > 0; --i)
			pmp[i][j] = pm - pd / (pmp[i + 1][j]) * pu;
	}

	//Upper and lower bounds 
	double lambda_L = (CorP == 'P') ? max(K - S[2 * NPaths], 0.0) : 0.0;
	double lambda_U = (CorP == 'C') ? max(S[0] - K, 0.0) : 0.0;

	// Work backwards and obtain matrix of Opt
	for (int j = NSteps; j > 0; --j)
	{
		pp[2 * NPaths - 1][j] = -pu * Opt[2 * NPaths - 2][j] - (pm - 2)*Opt[2 * NPaths - 1][j] - pd*Opt[2 * NPaths][j] + pd*lambda_L;

		for (int i = 2 * NPaths - 2; i > 0; --i)
		{
			pp[i][j] = -pu * Opt[i - 1][j] - (pm - 2)*Opt[i][j] - pd*Opt[i + 1][j] - pd / pmp[i + 1][j] * pp[i + 1][j];
		}

		for (int i = 0; i < 2 * NPaths + 1; ++i)
		{
			if (i == 0)
				C[i][j - 1] = (pp[i + 1][j] + pmp[i + 1][j] * lambda_U) / (pmp[i + 1][j] + pu);
			else if (i < 2 * NPaths)
				C[i][j - 1] = (pp[i][j] - pu*C[i - 1][j - 1]) / pmp[i][j];
			else
				C[i][j - 1] = C[i - 1][j - 1] - lambda_L;

			if (AorE == 'A')
				Opt[i][j - 1] = max((CorP == 'P') ? K - S[i] : S[i] - K, C[i][j - 1]);
			else
				Opt[i][j - 1] = C[i][j - 1];
		}

	}
	return Opt[NPaths][0];
}


//Explicit Finite Difference 
double Explicit_FD_S(double S0, double K, double T, double r, double sigma, char CorP, char AorE, int n = 250, int m = 40)
{
	// CorP: 'C' for Call option, 'P' for Put option
	// AorE: 'A' for American option, 'E' for European option

	double dS = double(2.0*K/m);
	int NPaths = m;				// Number of stock price paths
	double dt = double(T / n);				// Time increment
	int NSteps = n;		  // Number of time steps

	double pu ;         //Up probability
	double pm ;                  //Middle probability
	double pd ;     //Down probability

	  //Initialize stock price and option values.
	vector<double> S( NPaths + 1, 0.0);
	vector<vector<double>> Opt( NPaths + 1, vector<double>(NSteps + 1, 0.0));

	for (int j = 0; j <= NPaths; j++)
	{
		S[j] = double(j)*dS;
	}

	// Option price at maturity
	for (int i = 0; i < Opt.size(); i++)
		Opt[i][Opt[0].size() - 1] = ((CorP == 'P') ? max(K - S[i], 0.0) : max(S[i] - K, 0.0));

	// Work backwards and obtain matrix of Opt
	for (int j = NSteps - 1; j >= 0; j--)
	{
			for (int i = 1; i <  NPaths; ++i)
		{
			pu = 0.5*sigma*sigma*i*i*dt + 0.5*r*i*dt;
			pm = 1. - (sigma*sigma*i*i +r)*dt;
			pd = 0.5*sigma*sigma*i*i*dt - 0.5*r*i*dt;

			Opt[i][j] = (pu*Opt[i + 1][j + 1] + pm*Opt[i][j + 1] + pd*Opt[i - 1][j + 1]);
		}

			if (CorP == 'C') {
				//Lower boundary
				Opt[NPaths][j] = Opt[NPaths - 1][j] + (S[NPaths] - S[NPaths - 1]);
				//Upper boundary
				Opt[0][j] = Opt[1][j];
			}
			else if (CorP == 'P')
			{
				//Lower boundary
				Opt[NPaths][j] = Opt[NPaths - 1][j];
				//Upper boundary
				Opt[0][j] = Opt[1][j] - (S[1] - S[0]);
			}

		if (AorE == 'A') {
			for (int i = 0; i <=  NPaths; ++i)
				Opt[i][j] = ((CorP == 'P') ? max(K - S[i], Opt[i][j]) : max(S[i] - K, Opt[i][j]));
		}

	}
	int jstar = S0 / dS;
	double p = 0.;
	// run 2 point Lagrange polynomial interpolation
	p = p + (S0 - S[jstar + 1]) / (S[jstar] - S[jstar + 1])*Opt[jstar][0];
	p = p + (S0 - S[jstar]) / (S[jstar + 1] - S[jstar])*Opt[jstar+1][0];
	return p;
}


// Implicit Finite Difference
double  Implicit_FD_S(double S0, double K, double T, double r, double sigma, char CorP, char AorE, int n = 250, int m = 40)
{
	int NPaths = m;				// Number of stock price paths
	double dS = double(2.0*K / m);
	double dt = double(T / n);				// Time increment
	int NSteps = n;		  // Number of time steps

	double pu;      //Up probability
	double pm;      //Middle probability
	double pd;      //Down probability

					//Initialize stock price and option values.
	vector<double> S(NPaths + 1, 0.0);
	vector<vector<double>> Opt(NPaths + 1, vector<double>(NSteps + 1, 0.0));

	for (int j = 0; j <= NPaths; j++)
	{
		S[j] = double(j)*dS;
	}

	// Option price at maturity
	for (int i = 0; i < Opt.size(); i++)
		Opt[i][Opt[0].size() - 1] = ((CorP == 'P') ? max(K - S[i], 0.0) : max(S[i] - K, 0.0));

	//Initialize required matrices for solving the triangular system of equations.
	vector <vector<double>> pmp(NPaths, vector<double>(NSteps));
	vector <vector<double>> pp(NPaths, vector<double>(NSteps));
	vector<double> C(NPaths+1 );

	//Upper and lower bounds 
	double lambda_L = (CorP == 'P') ? max(K - S[NPaths], 0.0) : 0.0;
	double lambda_U = (CorP == 'C') ? max(S[0] - K, 0.0) : 0.0;

	//Solve the triangular system of equations.
	for (int j = NSteps - 1; j >= 0; --j)
	{
		for (int i = 0; i <= NPaths; ++i)
		{
			if (i == 0)
				C[i] = lambda_L;
			else if (i <  NPaths)
				C[i] = Opt[i][j + 1];
			else
				C[i] = lambda_U;
		}

		pm = 1. + (sigma*sigma*(NPaths - 1)*(NPaths - 1) + r)*dt;
		pu = -0.5*(sigma*sigma*(NPaths - 1)*(NPaths - 1) + r*(NPaths - 1))*dt;

		pmp[NPaths - 1][j] = pu + pm;
		pp[NPaths - 1][j] = C[NPaths - 1] + pu*lambda_U;
		for (int i = NPaths - 2; i > 0; --i)
		{
			pu = -0.5*(sigma*sigma*i*i + r*i)*dt;
			pm = 1. + (sigma*sigma*i*i + r)*dt;
			pd = -0.5*(sigma*sigma*i*i - r*i)*dt;

			pmp[i][j] = pm - pu*pd / (pmp[i + 1][j]);
			pp[i][j] = C[i] - pp[i + 1][j] * pu / pmp[i + 1][j];
		}

		pd = 0.;
		Opt[0][j] = (pp[1][j] + pmp[1][j] * lambda_L) / (pd + pmp[1][j]);
		Opt[1][j] = Opt[0][j] - lambda_L;

		for (int i = 2; i < NPaths; i++)
		{
			pd = -0.5*(sigma*sigma*i*i - r*i)*dt;
			Opt[i][j] = (pp[i][j] - pd*Opt[i - 1][j]) / pmp[i][j];
		}

		Opt[NPaths][j] = Opt[NPaths - 1][j] - lambda_U;

		if (AorE == 'A') {
			for (int i = 0; i <= NPaths; ++i)
				Opt[i][j] = ((CorP == 'P') ? max(K - S[i], Opt[i][j]) : max(S[i] - K, Opt[i][j]));
		}

	}
	int jstar = S0 / dS;
	double p = 0.;
	// run 2 point Lagrange polynomial interpolation
	p = p + (S0 - S[jstar + 1]) / (S[jstar] - S[jstar + 1])*Opt[jstar][0];
	p = p + (S0 - S[jstar]) / (S[jstar + 1] - S[jstar])*Opt[jstar + 1][0];
	return p;
}


//Crank Nicolson, BS PDE for S
double CrankNicolson_FD_S(double S0, double K, double T, double r, double sigma, char CorP, char AorE, int n = 250, int m = 40)

{
	double dS = double(2.0*K / m);
	int NPaths = m;				// Number of stock price paths
	int NSteps = n;		  // Number of time steps
	double dt = double (T/n);			// Time increment

	double pu ;      //Up probability
	double pm;      //Middle probability
	double pd ;      //Down probability
	 //Initialize stock price and option values.
	vector<double> S(NPaths + 1, 0.0);
	vector<vector<double>> Opt( NPaths + 1, vector<double>(NSteps + 1, 0.0));

	for (int j = 0; j <= NPaths; j++)
	{
		S[j] = double(j)*dS;
	}

	// Option price at maturity
	for (int i = 0; i < Opt.size(); i++)
		Opt[i][Opt[0].size() - 1] = ((CorP == 'P') ? max(K - S[i], 0.0) : max(S[i] - K, 0.0));

	//Initialize required matrices for solving the triangular system of equations.
	vector <vector<double>> pmp( NPaths , vector<double>(NSteps + 1));
	vector <vector<double>> pp(NPaths , vector<double>(NSteps + 1));
	vector <vector<double>> C( NPaths + 1, vector<double>(NSteps + 1));

	//Solve the triangular system of equations.
	for (int j = NSteps; j > 0; --j)
	{
		pm = 1. + 0.5*(sigma*sigma*(NPaths - 1)*(NPaths - 1) + r)*dt;
		pu = -0.25*(sigma*sigma*(NPaths - 1)*(NPaths - 1) + r*(NPaths - 1))*dt;
		pmp[ NPaths - 1][j] = pu + pm;

		for (int i = NPaths - 2; i > 0; --i)
		{
			pu = -0.25*(sigma*sigma*i*i + r*i)*dt;
			pm = 1. + 0.5*(sigma*sigma*i*i + r)*dt;
			pd = -0.25*(sigma*sigma*i*i - r*i)*dt;
			pmp[i][j] = pm - pu / (pmp[i + 1][j]) * pd;
		}
	}

	//Upper and lower bounds 
	double lambda_L = (CorP == 'P') ? max(K - S[ NPaths], 0.0) : 0.0;
	double lambda_U = (CorP == 'C') ? max(S[0] - K, 0.0) : 0.0;

	// Work backwards and obtain matrix of Opt
	for (int j = NSteps; j > 0; --j)
	{
		pu = -0.25*(sigma*sigma*(NPaths - 1 )*(NPaths - 1) + r*(NPaths - 1))*dt;
		pm = 1. + 0.5*(sigma*sigma*(NPaths - 1)*(NPaths - 1) + r)*dt;
		pd = -0.25*(sigma*sigma*(NPaths - 1)*(NPaths - 1) - r*(NPaths - 1))*dt;

		pp[NPaths - 1][j] = -pd * Opt[ NPaths - 2][j] - (pm - 2)*Opt[NPaths - 1][j] - pu*Opt[NPaths][j] + pu*lambda_U;

		for (int i =  NPaths - 2; i > 0; --i)
		{
			pu = -0.25*(sigma*sigma*i*i + r*i)*dt;
			pm = 1. + 0.5*(sigma*sigma*i*i + r)*dt;
			pd = -0.25*(sigma*sigma*i*i - r*i)*dt;
			pp[i][j] = -pd * Opt[i - 1][j] - (pm - 2)*Opt[i][j] - pu*Opt[i + 1][j] - pu / pmp[i + 1][j] * pp[i + 1][j];
		}

		for (int i = 0; i <  NPaths + 1; ++i)
		{
			pd = -0.25*(sigma*sigma*i*i - r*i)*dt;
			if (i == 0)
				C[i][j - 1] = (pp[i + 1][j] + pmp[i + 1][j] * lambda_L) / (pmp[i + 1][j] + pd);
			else if (i <  NPaths)
				C[i][j - 1] = (pp[i][j] - pd*C[i - 1][j - 1]) / pmp[i][j];
			else
				C[i][j - 1] = C[i - 1][j - 1] - lambda_U;

			if (AorE == 'A')
				Opt[i][j - 1] = max((CorP == 'P') ? K - S[i] : S[i] - K, C[i][j - 1]);
			else
				Opt[i][j - 1] = C[i][j - 1];
		}

	}
	int jstar = S0 / dS;
	double p = 0.;
	// run 2 point Lagrange polynomial interpolation
	p = p + (S0 - S[jstar + 1]) / (S[jstar] - S[jstar + 1])*Opt[jstar][0];
	p = p + (S0 - S[jstar]) / (S[jstar + 1] - S[jstar])*Opt[jstar + 1][0];
	return p;
}