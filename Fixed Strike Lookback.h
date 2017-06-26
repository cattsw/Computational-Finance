
// Fixed Strike Lookback Call and Put

#include <vector>
#include <math.h>
#include <random>
#include <iostream>
#include <algorithm>

using namespace std;

//Q1:  price of Lookback Call and Put by Monte Carlo simulation

double FSLookback(int seed, double S0, double K,  double r, double T, double sigma, int NSteps,int NPaths,char CorP)
{

	mt19937_64 rng;   // Random number generator and set the seed
	rng.seed(seed);				//set the seed	
	normal_distribution<double> N(0.0, 1.0);	// standard normal Z(0,1)
	double z; 
	double dt = T / double(NSteps); //time increasement
	vector< vector<double>> SPaths(NPaths, vector<double>(NSteps + 1));	
	vector<double> maxS(NPaths, 0.0);	 // Initialize maximum S for each path
	vector<double> minS(NPaths, 1e10);	 // Initialize minimum S for each path
	vector<double> Payoff(NPaths, 0.0);	 // Initialize  payoff
	double Opt=0.0;        //option price 

	for (int i = 0; i < NPaths; i++) {
		SPaths[i][0] = S0;
		for (int j = 1; j <= NSteps; j++)
		{
			z = N(rng);
			SPaths[i][j] = SPaths[i][j - 1] * exp((r - 0.5*sigma *sigma)*dt + sigma*sqrt(dt)*z);
			if (SPaths[i][j] > maxS[i]) maxS[i] = SPaths[i][j];
			if (SPaths[i][j] < minS[i]) minS[i] = SPaths[i][j];
		}
		}

		switch (CorP) {
		case 'C': {
			for (int i = 0; i < NPaths; i++) {
				Payoff[i] = max(maxS[i] - K, 0.0);
				Opt = Opt + exp(-r*T)*Payoff[i];
			}
		}
			 break;
		case 'P': {
			for (int i = 0; i < NPaths; i++) {
				Payoff[i] = max(K - minS[i], 0.0);
				Opt = Opt + exp(-r*T)*Payoff[i];
			}
		}
			 break;
		default:   cout << "Error, choose C or P \n"; break;   // wrong input
		
	   }
		
	return Opt / NPaths;
}


