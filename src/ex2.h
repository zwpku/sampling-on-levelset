#include "ranlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <string.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <assert.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std ;

extern "C" {
  extern void dgeqrf_(int *, int * , double *, int * , double *, double *, int * , int *) ;

  extern void dorgqr_(int *, int * , int *, double *, int * , double *, double *, int * , int *) ;

  extern void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *) ;

  extern void dgetrf_(int *, int *, double *, int *, int *, int *) ;
}

const double pi = atan(1) * 4 ;

extern int n, N, d, k, n_bins , newton_max_step, verbose_flag , output_every_step, *ipiv, max_lag ;

extern int maximal_try_number, mcmc_flag ;

extern double trace_b, h, h_mcmc, eps_tol, *mat_x, determinant_tol, mean_xi_distance, bin_width ;

extern ofstream log_file ;

// for SDE-based simulation
extern int maximal_try_number ;
// step-size in ODE solver
extern double dt ;

extern vector<vector<double> > grad_vec;

extern vector<int> counter_of_each_bin ;

extern vector<double> trace_series ;

void xi( vector<double> & x, vector<double> & ret ) ;

double trace(vector<double> & x) ;

int check_determinant(vector<double> & x) ;

void grad_xi(vector<double> & x, vector<vector<double> > & grad) ;

void grad_xi(vector<double> & x, vector<vector<double> > & grad) ;

double vec_dot(vector<double> & v1, vector<double> & v2) ;

double distance_to_levelset(vector<double> & x) ;

void init_rand_generator() ;

void analysis_data_and_output() ;
