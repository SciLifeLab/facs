extern double get_mu (long long num_hit, double prob);
extern double get_sigma (long long num_hit, double prob);
double erf (double x);
double pdf (double x, double mu, double sigma);
extern double cdf (double x, double mu, double sigma);
extern double prob_suggestion (int k_mer);
