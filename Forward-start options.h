
// Forward start option pricing

using namespace std;


// Monte Carlo simulation for pricing European forward start options

double ForSt_EuPut(double S0,  double r,  double sigma,double t1, double T,int NSteps, int n )
{
	// Stock paths up to T-t1
	vector<vector<double>> S= Stock_Path(54321, 1.0, r, sigma, T-t1, NSteps, n);
	int NPaths = 2 * n;  // antithetic 
	vector<double> Payoff(NPaths, 0.0);	// ATM option price at future time t1.
	double EuP=0.0;        //at the money European put price with s=1
	double K=1.0;      //strike price

	for (int i = 0; i <  NPaths; i++) {
		
		Payoff[i] = max(K-S[i][NSteps] , 0.0);
		EuP += exp(-r*(T-t1))*Payoff[i];
	}

	return S0*EuP / NPaths;
}


//   Monte Carlo simulation  for pricing American forward start options

double ForSt_AmPut(double S0, double r, double sigma, double t1, double T, int NSteps, int n)
{
	//Stock paths before start time
	vector<vector<double>> S = Stock_Path(123456, 1.0, r, sigma, T-t1, NSteps, n);
	
	int NPaths = 2 * n;  // antithetic 
	double K=1.0;      //strike price		   
	double AmPut = LSMC(S, K, r, T - t1, NSteps, NPaths, 'L', 3, 0); //  ATM American option price at future time t1.
	
	return S0*AmPut;
}