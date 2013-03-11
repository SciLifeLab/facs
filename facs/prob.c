//#include <iostream>
//#include <cmath>
//using namespace std;

#include <stdio.h>
#include <math.h>
#include "prob.h"

// Returns the erf() of a value (not super precice, but ok)
//double erf(double x);
//double pdf(double x, double mu, double sigma);
//extern double cdf(double x, double mu, double sigma);
//int main();


double erf(double x)
{  
 double y = 1.0 / ( 1.0 + 0.3275911 * x);   
 return 1 - (((((
        + 1.061405429  * y
        - 1.453152027) * y
        + 1.421413741) * y
        - 0.284496736) * y 
        + 0.254829592) * y) 
        * exp (-x * x);      
}

// Returns the probability of x, given the distribution described by mu and sigma.
double pdf(double x, double mu, double sigma)
{
  //Constants
  static const double pi = 3.14159265; 
  return exp( -1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * pi));
}

// Returns the probability of [-inf,x] of a gaussian distribution
double cdf(double x, double mu, double sigma)
{
	return 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2.))));
}

/*
int main()
{
	
	double x, mu, sigma;
	cout << "x, mu, sigma: ";
	cin >> x >> mu >> sigma;

	cout << "PDF of x is: " << pdf(x,mu,sigma) << endl;
	cout << "CDF of x is: " << cdf(x,mu,sigma) << endl;

        x = 1;
        mu = 10;
        sigma = 2.5;
        printf ("CDF->%e\n",cdf(x,mu,sigma));

	return 0;
}
*/
