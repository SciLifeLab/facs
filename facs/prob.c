#include <stdio.h>
#include <math.h>
#include "prob.h"  
//double erf(double x);
//double pdf(double x, double mu, double sigma);
//extern double cdf(double x, double mu, double sigma);
double erf (double x) 
{
  double y = 1.0 / (1.0 + 0.3275911*x);
  return 1 - (((((+1.061405429*y-1.453152027)*y+1.421413741)*y-0.284496736)*y+0.254829592)*y)*exp(-x*x);
}
// Returns the probability of x, given the distribution described by mu and sigma.
double pdf (double x, double mu, double sigma) 
{ 
  //Constants
  static const double pi = 3.14159265;
  return exp (-1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma *sqrt (2 *pi));
}
// Returns the probability of [-inf,x] of a gaussian distribution
double cdf (double x, double mu, double sigma) 
{
  return 0.5 * (1 + erf ((x - mu) / (sigma * sqrt (2.))));
}

double get_mu (long long num_hit, double prob)
{
  return ((double) num_hit) * prob;
}

double get_sigma (long long num_hit, double prob)
{
  return (double) num_hit *prob * (1 - prob);
}
/*
int main()
{	
        printf ("CDF->%e\n",cdf(18783005,get_mu(18783005000,0.026),get_sigma(18783000005,0.026)));
	return 0;
} 
*/
