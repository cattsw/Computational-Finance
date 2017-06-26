//CMF Project2 by Shunwen Tan

//function of generating random numbers

#include <vector>
#include <limits>
#include <cstddef>
#include <math.h>
using namespace std;

// Independent normally distributed random numbers 
double Random_Normal(double mu = 0.0, double sigma = 1.0)
{
	const double epsilon = (std::numeric_limits<double>::min)();
	const double pi = 3.14159265358979323846;
	double u1, u2,z;

	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	} while (u1 <= epsilon);   // control for log

	z = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2);
	return mu+z*sigma;
}

// Normally distributed random numbers with correlation a
vector<double> Random_CorNormal(double mu1 = 0.0, double sigma1 = 1.0, double mu2 = 0.0, double sigma2 = 1.0, double a = -0.7)
{

	double x0, x1;
	double z0, z1;
	vector<double>  y;

	z0 = Random_Normal();
	z1 = Random_Normal();
	x0 = mu1 + z0*sigma1;
	y.push_back(x0);
	x1 = mu2 + sigma2 *(a*z0 + sqrt(1.0 - a*a)*z1);
	y.push_back(x1);

	return y;
}