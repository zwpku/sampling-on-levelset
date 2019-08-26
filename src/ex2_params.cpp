#include "ex2.h"

int n, N, d, k ;
int n_bins, output_every_step, maximal_try_number ;

int newton_max_step, max_lag ;

double h, dt, eps_tol, h_mcmc ;

double bin_width, trace_b ;

double mean_xi_distance ;

double mean_trace ;

int verbose_flag, mcmc_flag ;
ofstream log_file ;

vector<int> counter_of_each_bin ;

vector<double> trace_series ;

vector<vector<double> > grad_vec;

double * mat_x ; 
int *ipiv ;

double determinant_tol ;

