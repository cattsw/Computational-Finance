
// Apply jump diffusion to pricing of default options

#include <vector>
#include <math.h>
#include <random>
#include <iostream>
#include <algorithm>
using namespace std;

//default option 
vector<double> Proj6_2func( double lambda1, double lambda2, double T)
{
	double V0 = 20000.0;   //initial value of collateral
	double L0 = 22000.0;  //initial outstanding balance of loan
	// K  = strike price
	double mu = -0.1;  // drift of collateral
	double sigma = 0.2; // volatillity of collateral		
	double dt = 0.01;// Time increment
	int NSteps = int(T/dt); // number of time steps
	int NPaths = 10000; // number of simulations

	double gamma = -0.4; //jump parameter
	
	//loan parameters
	double r0 = 0.02;
	double delta = 0.25;
	double R = r0 + delta*lambda2;//APR of the loan
	double r = R / 12.0;  //contract rate of loan
	double n = T * 12.0;
	double PMT = L0*r / (1.0 - 1.0 / pow(1.0 + r, n));
	double a = PMT / r;
	double b = PMT / (r*pow(1.0 + r, n));
	double c = 1.0 + r;

	//default parameters
	double alpha = 0.7;
	double epsilon = 0.95; //recovery rate of the collateral
	double beta = (epsilon - alpha) / T;
	


	mt19937_64 rng;   // Random number generator and set the seed
	rng.seed(1234);		//set the seed	

	if (lambda1 == 0.0) lambda1 = 0.0000000000001;
	if (lambda2 == 0.0) lambda2 = 0.0000000000001;
	poisson_distribution<int> P1(lambda1*dt); 	// Define the poisson distributions for jump
	poisson_distribution<int> P2(lambda2*dt); 	// Define the poisson distributions for stopping time
	normal_distribution<double> NS(0.0, 1.0);  	// Define the standard normal distributions

	// Initialize the option payoff ,number of defaluts, default time
	double Payoff = 0.0;
	int DefN = 0;
	double DefT=0.0;

	// Initialize the collateral and loan paths 
	vector<vector<double> > V(NPaths, vector<double>(NSteps + 1));
	vector<vector<double> > L(NPaths, vector<double>(NSteps + 1));
	vector<vector<double> > q(NPaths, vector<double>(NSteps + 1));

	for (int s = 0; s<NPaths; s++) {
		V[s][0] = V0;
		L[s][0] = L0;
		for (int t = 1; t<= NSteps; t++) {
			t = double(t);
	/*		double J = 0.0;
			if (lambda1 >0.0000000000001) {
				int Nt1 = P1(rng);
				if (Nt1 > 0)
					for (int i = 0; i<Nt1; i++)
						J += gamma;
			} 
			double z = NS(rng);
			V[s][t] = V[s ][t-1] * exp((mu - 0.5*sigma*sigma)*dt + sigma*sqrt(dt)*z + J);	
	*/		
			double J = 1.0;
			if (lambda1 >0.0000000000001) {
				int Nt1 = P1(rng);
				if (Nt1 > 0)
					for (int i = 0; i<Nt1; i++)
						J *= (1+gamma);
			}
			double z = NS(rng);
			V[s][t] = V[s ][t-1] * exp((mu - 0.5*sigma*sigma)*dt + sigma*sqrt(dt)*z )*J;
			
			L[s][t] = a - b*pow(c, 12.0*t*dt);  //outstanding balance of loan
			q[s][t] = alpha + beta*t*dt;
	
			 if (V[s][t] <= q[s][t] * L[s][t]) {
				Payoff += exp(-r0*t*dt)*max(L[s][t] - epsilon*V[s][t],0.0);
				DefN += 1;
				DefT += t*dt;
		//		cout << s << "," << t << endl;
				break;
			}
			 else
			 {
				 int Nt2 = 0;
				 if (lambda2 > 0.0000000000001)  Nt2 = P2(rng);
				 if (Nt2 > 0) {
					 Payoff += exp(-r0*t*dt)*abs(L[s][t] - epsilon*V[s][t]);
					 DefN += 1;
					 DefT += t*dt;
					 //		cout << s << "," << t << endl;
					 break;
				 }
			 }
		//	cout << s << "," << "no exercise" << endl;
		}
	}
	double D= Payoff / NPaths;
	double DProb=double(DefN)/ NPaths;
	double EET=1000000;   
	if (DefN > 0) EET = DefT / DefN;
	
	vector<double> output(3);
	output[0] = D;
	output[1] = DProb;
	output[2] = EET;
	return output;
}


