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

const double pi = atan(1) * 4 ;

int n ;
int tot_step ;
int n_bins ;

double h, beta , noise_coeff ; 

double eps_tol, reverse_tol, newton_grad_tol ;
double mean_xi_distance ;
double pot_coeff ; 

int newton_max_step ;

int forward_newton_counter, metropolis_counter, backward_newton_counter, reverse_check_counter ;

/*
 * c2=c^2, c4=c^4
 */
double c, c2, c4 ;
double bin_width ;

vector<double> state, tangent_v, tmp_state, grad_vec, counter_of_each_bin ;

// the reaction coordinate function 
double xi(vector<double> & x)
{
  return 0.5 * (x[0] * x[0] / c2 + x[1] * x[1] - 1.0) ;
}

// gradient of the reaction coordinate function 
void grad_xi(vector<double> & x, vector<double> & grad)
{
  grad[0] = x[0] / c2 ;
  grad[1] = x[1] ;
}

double U(vector<double> & x)
{
  double tmp ;
  tmp = x[0] - c ;
//  tmp = x[1] - 1.0 ;
  return  pot_coeff * tmp * tmp * 0.5 ;
}

double grad_U(vector<double> & x, vector<double> & grad)
{
  double tmp ;
  tmp = x[0] - c;
//  tmp = x[1] - 1.0 ;
  grad[0] = pot_coeff * tmp ; grad[1] = 0.0 ;
}

double vec_dot(vector<double> & v1, vector<double> & v2) 
{
  return v1[0] * v2[0] + v1[1] * v2[1] ;
}

double len_grad_xi(vector<double> & x)
{
  vector<double> tmp_vec ;
  tmp_vec.resize(2) ;
  grad_xi(x, tmp_vec) ;
 // return 1.0 ;
  return sqrt(vec_dot(tmp_vec, tmp_vec)) ;
}

/*
 * Compute the projection by Newton-method.
 *
 */
int newton_projection(vector<double> & state, vector<double>& vec, vector<double> & result)
{
  double eps, a, tmp ;
  int step ;
  vector<double> k1, k2 ;

  k1.resize(2) ;  k2.resize(2) ;  

  a = 0 ; 
  step = 0 ; 
  k1 = state ;
  eps = fabs(xi(k1)) ;

  while (eps > eps_tol)  
  {
    grad_xi(k1, k2) ;
    tmp = vec_dot(k2, vec) ;

    // unsuccessful, if gradient is too small
    if (fabs(tmp) < newton_grad_tol) 
    {
      tot_step += step ;
      return 0 ;
    }

    a = a - xi(k1) / tmp ;

    step ++ ;
    // unsuccessful, if maximal steps are reached
    if (step > newton_max_step) // not success
    {
      tot_step += step ;
      return 0 ;
    }

    k1[0] = state[0] + a * vec[0] ;
    k1[1] = state[1] + a * vec[1] ;
    eps = fabs(xi(k1)) ;
  }

  tot_step += step ;
  result = k1 ;

//  printf("step=%d\teps=%.4e\t %.4e\ta=%.4f\n", step, fabs(xi(state)), eps, a) ;
  // successful 
  return 1; 
}

double dist_xy(vector<double> & x, vector<double> & y)
{
  return sqrt( (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) ) ;
}

/*
 * update state by MCMC
 *
 */
void update_state(vector<double> & state )
{
  double r, norm , r1, accept_prob , tmp ;
  int flag ;
  vector<double> y_state, yy_state ;

  y_state.resize(2) ;
  yy_state.resize(2) ;
  
  /* 
   * normal direction 
   */
  grad_xi(state, grad_vec) ;

  norm = sqrt( vec_dot(grad_vec, grad_vec) ) ;

  // unit tangent vector 
  tangent_v[0] = - grad_vec[1] / norm ;
  tangent_v[1] = grad_vec[0] / norm ;

  // generate standard Guassian variable 
  r = gennor(0, 1) ;

  // a move along tangent vector v
  for (int i = 0 ; i < 2 ; i ++)
  {
    tmp_state[i] = state[i] + noise_coeff * r * tangent_v[i] ;
  }

  mean_xi_distance += fabs(xi(tmp_state)) ;

  /* 
   * Step 2: projection by Newton-method.
   */
  flag = newton_projection(tmp_state, grad_vec, y_state) ;

  if (flag == 0) {
    forward_newton_counter ++ ;
    return ;
  }

  /* 
   * normal direction at y 
   */
  grad_xi(y_state, grad_vec) ;
  norm = sqrt( vec_dot(grad_vec, grad_vec) ) ;
  // unit tangent vector 
  tangent_v[0] = - grad_vec[1] / norm ;
  tangent_v[1] = grad_vec[0] / norm ;
  r1 = (vec_dot(tangent_v, state) - vec_dot(tangent_v, y_state)) / noise_coeff ;
  accept_prob = exp(-beta * (U(y_state) - U(state))) * len_grad_xi(state) / len_grad_xi(y_state) * exp((r * r - r1 * r1) * 0.5) ;
  if (accept_prob > 1.0) accept_prob = 1.0 ;
  tmp = ranf() ;
  if (tmp > accept_prob) 
  {
    metropolis_counter ++ ;
    return ;
  }

  // a move along tangent vector v
  for (int i = 0 ; i < 2 ; i ++)
  {
    tmp_state[i] = y_state[i] + noise_coeff * r1 * tangent_v[i] ;
  }

  flag = newton_projection(tmp_state, grad_vec, yy_state) ;

  if (flag == 0)
  {
    backward_newton_counter ++ ;
    return ;
  }

  if (dist_xy(yy_state, state) > reverse_tol) 
  {
    reverse_check_counter ++ ;
    return ;
  }

  state = y_state ;
}

/* 
 * Initialize the seed of random numbers depending on the cpu time 
 */
void init_rand_generator()
{
  long is1, is2 ;
  long current_time ;

  char phrase[100] ;

  current_time = time(NULL) ;

  sprintf( phrase, "%ld", current_time ) ;

  phrtsd(phrase, &is1, &is2) ;
  setall(is1, is2) ;
}

int main ( int argc, char * argv[] ) 
{
  char buf[50] ;
  ofstream out_file ;
  int idx ;
  int output_every_step ;
  double angle , T ;

  clock_t start , end ;

  n = 5000000 ;
  h = 5.00 ;
  T = n * h ;
  pot_coeff = 10.0 ;

  // compute the total steps
  output_every_step = 1000 ;

  beta = 1.0 ;

  c = 3.0 ;
  c2 = c*c ;
  c4 = c2 * c2 ;

//  eps_tol = 1e-7 ;
  eps_tol = 1e-12 ;
  newton_grad_tol = 1e-7 ; 
  newton_max_step = 10 ;
  reverse_tol = 1e-10 ;

  tot_step = 0 ; 
  forward_newton_counter = 0 ;
  backward_newton_counter = 0 ;
  metropolis_counter = 0 ;
  reverse_check_counter = 0 ;

  mean_xi_distance = 0 ;

  n_bins = 50 ;

  // divied [0, 2pi] to n_bins with equal width
  bin_width = 2 * pi / n_bins ;

  counter_of_each_bin.resize(n_bins,0) ;

  init_rand_generator();
  noise_coeff = sqrt(2.0 / beta * h) ;

  //initial state
  state.resize(2) ;
  tangent_v.resize(2) ;
  grad_vec.resize(2) ;
  tmp_state.resize(2) ;
  state[0] = 0 ; state[1] = 1.0 ;

  start = clock() ;

  printf("Total time = %.2f\nh=%.2e\nn=%.1e\nNo. of output states=%d\n", T, h, n *1.0, n / output_every_step) ;

  sprintf(buf, "../data/traj_mcmc.txt" ) ;

  out_file.open(buf) ;

  out_file << n / output_every_step << endl ;

  angle = 0 ;
  for (int i = 0 ; i < n ; i ++)
  {
    if (i % output_every_step == 0)
      out_file << state[0] << ' ' << state[1] << endl ;

    update_state(state) ;

    //compute histgram during the simulation
    angle = atan2(state[1], state[0] / c) ;
    // change angle to [0, 2*pi]
    if (angle < 0) angle = 2 * pi + angle ;
    idx = int (angle / bin_width) ;
    counter_of_each_bin[idx] ++ ;
  }

  int tot_rej ;
  out_file.close() ;
  printf("\naverage iteration steps = %.2f\n", tot_step * 1.0 / n) ;
  printf("\naverage xi distance = %.3e\n", mean_xi_distance * 1.0 / n) ;
  tot_rej = forward_newton_counter + backward_newton_counter + reverse_check_counter + metropolis_counter ;
  printf("\nRejection rate: Forward\tReverse\tReversibility\tMetrolis\tTotal\n") ;
  printf("\t\t %.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", forward_newton_counter * 1.0 / n, backward_newton_counter * 1.0 / n, reverse_check_counter *1.0/n, metropolis_counter * 1.0 / n, tot_rej * 1.0 / n) ;

  sprintf(buf, "../data/counter_mcmc.txt") ;
  out_file.open(buf) ;

  out_file << n << ' ' << n_bins << endl ;

  for (int i = 0 ; i < n_bins ; i++)
  {
    out_file << counter_of_each_bin[i] << ' ';
  }
  out_file << endl ;
  out_file.close() ;

  end = clock() ;
  printf("\n\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;

  return 0; 
}
