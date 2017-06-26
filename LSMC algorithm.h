// Stock paths, LSM algorithm

#include "StdAfx.h"
#include <iostream>
#include <math.h>
#include <random>
#include <algorithm>
#include <vector>
#include <set>
#include <iterator>
#include <numeric>
#include <limits>
#include "MatrixOperations.h"
using namespace std;

// Generate 2*n antithetic Stock price paths 
vector<vector<double>> Stock_Path(int seed, double S0 , double r , double sigma, double T , int NSteps, int NPaths)
{
	srand(seed);
	double u1, u2;
	double z1, z2; //antithetic normal random 
	const double pi = 3.14159265358979323846;
	vector< vector<double>> SPaths(2*NPaths, vector<double>(NSteps + 1));
	double dt = T / double(NSteps); //step size
	// Floor u1 to avoid errors with log function
	const double epsilon = (std::numeric_limits<double>::min)();

	for (int i = 0; i < NPaths; i++)
	{
		SPaths[2*i][0] = S0;
		SPaths[2*i+1][0] = S0;
		for (int j = 1; j <= NSteps; j++)
		{
			// Independent uniform random variables
			do
			{
				u1 = rand() * (1.0 / RAND_MAX);
				u2 = rand() * (1.0 / RAND_MAX);
			} while (u1 <= epsilon);
			// Floor u1 to avoid errors with log function
			// Z ~ N(0,1) by Box-Muller transformation
			z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2);
			z2 = -z1;
			SPaths[2 * i][j] = SPaths[2 * i][j - 1] * exp((r - 0.5*sigma *sigma)*dt + sigma*sqrt(dt)*z1);
			SPaths[2 * i + 1][j] = SPaths[2 * i + 1][j - 1] * exp((r - 0.5*sigma *sigma)*dt + sigma*sqrt(dt)*z2);
		}
	}
	return SPaths;
}

// Basis functions with different orders
// Hermite polynomial
double Hermite(double x, int k) 
{
	double y=1.0;
	 if (k == 2) { y = 2.0*x; }
	else if (k==3) { y = 4.0*x*x-2.0; }
	else if (k == 4) { y = 8.0*x*x*x - 12.0 * x; }
	else if (k >= 5) {
		y = 2.0*x*Hermite(x, k - 1) - 2.0* double(k - 1)*Hermite(x, k - 2);}
	return y;
}

double Weight_Hermite(double x, int k)
{
	double y = 1.0;
	if (k >= 1)
		y = exp(-x*x / 2.0)*Hermite(x, k);
	return y;
}

// Laguerre polynomial 
double Laguerre(double x, int k)
{
	double y = 1.0;
	if (k == 2) { y = 1.0-x; }
	if (k == 3) { y = 1.0 - 2.0*x + x*x/2.0; }
	if (k==4) { y = 1.0 - 3.0*x + 3.0*x*x / 2.0-x*x*x/6.0; }
	else if (k >= 4) {
		y = 1.0/(k-1)*((2.0*k-3.0-x)*Laguerre(x,k-1)-double(k-2)*Laguerre(x, k - 2));
	}
	return y;
}

double Weight_Laguerre(double x, int k)
{
	double y = 1.0;
	if (k >= 1)
	y=exp(-x / 2.0)*Laguerre(x, k );
	return y;
}

// Monomials 
double Monomial(double x, int k)
{
	double y = pow(x,k-1);
	return y;
}

//LSMC algorithm for American Put option 
double LSMC(vector<vector<double>> S, double K, double r,  double T, int NSteps, int NPaths, char Fbasis,int fk,int W)
{
	// Fbasis = different types of polynormials: L,H,M
	// fk= number of basis functions
	// if W==1 weighted  polynormials,if W==0 unweighted

	double dt = T / double(NSteps); // Time increment
	vector<vector<double>> CF(NPaths, vector<double> (NSteps+1)); 	// Initialize the Cash Flows.

	// Last cash flows
	for (int i = 0; i < NPaths; i++)
			CF[i][NSteps] = max(K - S[i][NSteps], 0.0);

	// Work backwards through the stock prices until time t=1
	for (int t = NSteps - 1; t >= 1; t--)
	{
		// Indices for stock paths in-the-money at time t
		vector<int> Index(NPaths);
		for (int s = 0; s < NPaths; s++)
		{
			Index[s] = 0;
			if (S[s][t] < K) Index[s] = 1;
		}

		// Stock paths in-the-money at time t
		int NI = 0;
		vector<double> X;
		vector<int> X_I;
		for (int s = 0; s < NPaths; s++)
			if (Index[s] = 1) {
				X.push_back(S[s][t]);
				X_I.push_back(s);
				NI += 1;
			}

		// Discounted cash flows at time t+1
		vector<double> Y(NI);
		for (int s = 0; s < NI ; s++)
			Y[s] = CF[X_I[s]][t + 1] * exp(-r*dt);

		vector<vector<double> > Z(NI, vector<double>(fk));

		if (W == 0)
		{
			// Matrix for predicting regression base on different basis functions
			switch (Fbasis) {
			case 'H': {
				for (int s = 0; s <= NI - 1; s++) {
					for (int k = 1; k <= fk; k++) {
						Z[s][k - 1] = Hermite(X[s], k);
					}
				}
			}
					  break;
			case 'L': {
				for (int s = 0; s <= NI - 1; s++) {
					for (int k = 1; k <= fk ; k++) {
						Z[s][k-1] = Laguerre(X[s] , k);
					}
				}
			}
					  break;
			case 'M': {
				for (int s = 0; s <= NI - 1; s++) {
					for (int k = 1; k <= fk; k++) {
						Z[s][k - 1] = Monomial(X[s], k);
					}
				}
			}
					  break;
			default:            // wrong input
				cout << "Error, choose H/L/M \n";
				break;
			}
		}

		else if (W == 1)
		{
			// Matrix for predicting regression base on different basis functions
			switch (Fbasis) {
			case 'H': {
				for (int s = 0; s <= NI - 1; s++) {
				//	Y[s] = Y[s] / K;
					for (int k = 0; k <= fk-1; k++) {
						Z[s][k ] = Weight_Hermite(X[s] / K, k);//weighted Hermite, renormlize
					}
				}
			}
					  break;
			case 'L': {
				for (int s = 0; s <= NI - 1; s++) {
				//	Y[s] = Y[s] / K;
					for (int k = 0; k <= fk - 1; k++) {
						Z[s][k] = Weight_Laguerre(X[s] / K, k);//weighted Laguerre, renormlize
					}
				}
			}
					  break;
			default:            // wrong input
				cout << "Error, choose H or L \n";
				break;
			}
		
		}
		else printf("W should be set to 1(weighted) or 0(unweighted) \n");

		// Regression parameters and predicted cash flows
		vector<double> beta=Beta(Z, Y); 
		vector<double> PredCF = Mult_Vec(Z, beta);

		// Indices for stock paths that  exercise is optimal
		// E_I[s] contains the path number
		vector<int> E;
		vector<int> E_I;
		int NE = 0;
		for (int s = 0; s <= NI - 1; s++)
			if (( (K - X[s]>0)  &  (K - X[s]> PredCF[s])) )
			{
				E_I.push_back(s);
				E.push_back(X_I[s]);
				NE += 1;
			}

		// Indices for stock paths that continuation is optimal
		vector<int> All(NPaths);
		for (int k = 0; k < NPaths; k++)
			All[k] = k;

		// C contains indices for continuation paths: C = All - E 
		set<int> CSet;
		set_difference(All.begin(), All.end(), E.begin(), E.end(),inserter(CSet, CSet.end()));
		int NC = CSet.size();
		// Convert to a vector
		vector<int> C(CSet.begin(), CSet.end());

		// Replace cash flows with exercise value where exercise is optimal
		for (int s = 0; s < NE ; s++)
				CF[E[s]][t] = max(K - X[E_I[s]], 0.0);	
		for (int s = 0; s < NC ; s++)
			CF[C[s]][t] = exp(-r*dt)*CF[C[s]][t + 1];
		
	}

	//  American price
	double AmPut = 0.0;

	for (int s = 0; s < NPaths; s++) 
		AmPut += CF[s][1];

	return exp(-r*dt)*AmPut/NPaths;
}
