//MBS pricing with different CPR models

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <algorithm> 
#include "CIR model.h"
#include <iostream>

using namespace std;

//Numerix prepayment model for CPR

double CPR_Numerix(int t, double PVt, double PV0,double WAC, double rt , double rbar = 0.08, double kappa = 0.6, double sigma = 0.12)
{
	double r = WAC / 12.;
	double R = 12. * r;
	double RI; //refinance incentives
	double BU; // burnout
	double SG; //seasoning
	double SY; //seasonality
	vector<double> SYt = { 0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.1,1.18,1.22,1.23,0.98 };

	//RI
	double rT = -0.1*log(DBond_CIR(double(t - 1) / 12, double(t - 1) / 12 + 10., rt, rbar, kappa, sigma));

	RI = 0.28 + 0.14*atan(-8.57 + 430.*(R - rT));
	BU = 0.3 + 0.7*PVt / PV0;
	SG = min(1., t / 30.);
	SY = SYt[(t - 1) % 12];
	double CRP = RI*BU*SG*SY;
	return CRP;
}


//PSA model for CPR
double CPR_PSA(int t)
{
	vector<double> CPR(31, 0.);
	for (int i = 1; i <= 30; i++)
			CPR[i] = CPR[i-1] + 0.002;
	if (t <= 30) return CPR[t];
	else return CPR[30];
}


// price of MBS, IOs, POs
vector<double> MBS(char model, double kappa, double rbar, double x = 0.,
	double T = 30., double PV0 = 100000, double WAC = .08, double r0 = 0.078, double sigma = 0.12)
{
	// model= 'N': numerix-prepayment model, model= 'P': PSA model
	// x: OAS spread

	double r = WAC / 12.;
	double R = 12. * r;
	int n = 12 * T;// number of CFs
	double dt = 1. / 360;
	vector<double> CF(n + 1, 0.);
	vector<double> TPP(n + 1, 0.);
	vector<double> IP(n + 1, 0.);
	vector<double> PV(n + 1, 0.);
	vector<double> CPR(n + 1, 0.);
	double MBS = 0., IOs=0., POs=0.;
	PV[0] = PV0;

	for (int t = 1; t <= n; t++)
		{
		t = double(t);
		if (model == 'N') {
			//double rt = r0*exp(-kappa*(t / 12.)) + rbar*(1 - exp(-kappa*(t / 12.)))+x;
			double rt = r0*exp(-kappa*(t / 12.)) + rbar*(1 - exp(-kappa*(t / 12.))) ;
			CPR[t] = CPR_Numerix(t, PV[t - 1], PV0, WAC, rt, rbar, kappa, sigma);
		}
			else if (model == 'P')
				CPR[t] = CPR_PSA(t);

			IP[t] = PV[t - 1] * r;
			CF[t] = PV[t - 1] * r / (1. - pow(1. + r, -n + t - 1)) + (PV[t - 1] - PV[t - 1] * r*(1. / (1. - pow(1. + r, -n + t - 1)) - 1.))*(1. - pow(1. - CPR[t], 1. / 12));
			TPP[t] = CF[t] - IP[t];
			PV[t] = PV[t - 1] - TPP[t];

			double d_t =-log(DBond_CIR(0., t/12., r0, rbar, kappa, sigma))+x*t/12.;
		
		MBS = MBS + exp(-d_t)*CF[t];
		IOs = IOs + exp(-d_t)*IP[t];
		POs = POs + exp(-d_t)*TPP[t];
	}
	vector<double> output = { MBS , IOs , POs  };
	return output;
}

// Bisection Algorithm for root search x
double Root_x(double a, double b,double P)
{
	const int n = 1000000;  //number of max iterations
	double epsilon = 0.00001;   // Tolerance
	double mid = 0.5 * (a + b);
	double mid_P = 0., l_P = 0., h_P = 0.;
		mid_P = MBS('N',0.6,0.08, mid)[0] ;
		l_P = MBS('N', 0.6, 0.08, a)[0];
		h_P = MBS('N', 0.6, 0.08, b)[0];
	
	// make sure choose the appropriate interval
	double  l = P - l_P;
	double h = P - h_P;
	if (l*h > 0)
	{
		printf("wrong interval\n");
		return -1;
	}

	// While the difference between midP and the market price is greater than epsilon, keep subdividing the interval 
	int i = 0;
	do {
		if (mid_P < P) {
			a = mid;
		}
		else
		{
			b = mid;
		}
		mid = 0.5 * (a + b);
		
		mid_P = MBS('N', 0.6, 0.08, mid)[0];
		
		i = i + 1;

	} while ((abs(mid_P - P) > epsilon) & (i <= n));

	return mid;
}

