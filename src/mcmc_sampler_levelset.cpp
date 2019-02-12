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
}

const double pi = atan(1) * 4 ;

int n ;
int N, d, k ;
int tot_newton_step ;
int n_bins ;

double h, beta , noise_coeff ; 

double eps_tol, reverse_tol, newton_grad_tol ;
double mean_xi_distance ;

int newton_max_step ;

int forward_newton_counter, metropolis_counter, backward_newton_counter, reverse_check_counter ;

double bin_width ;

vector<double> state, tangent_v, tmp_state, grad_vec, counter_of_each_bin ;

// arrays for QR decomposition
double * qr_tau, * qr_work , * qr_mat ;

double * mat_a, * linear_sol; 
int *ipiv ;

// 
// the vector-valued reaction coordinate function: N + N * (N-1) / 2
//
void xi( vector<double> & x, vector<double> & ret )
{
  int idx, i0, j0 ;
  // first N component 
  for (int i = 0 ; i < N; i ++)
  {
    ret[i] = 0 ;
    i0 = i * N ;
    for (int j = 0 ; j < N; j ++)
      ret[i] += x[i0+j] * x[i0+j] ;
    ret[i] = 0.5 * (ret[i] - 1) ;
  }

  // the remaining N * (N-1) / 2 component
  idx = N ;
  for (int i = 0 ; i < N ; i ++)
    for (int j = i+1 ; j < N; j ++)
    {
      ret[idx] = 0 ;
      i0 = i * N; j0 = j * N ;
      for (int j1 = 0 ; j1 < N; j1 ++)
	ret[idx] += x[i0 + j1] * x[j0 + j1] ;
      idx ++ ;
    }
}

// gradient of the reaction coordinate functions
//
// grad is a matrix of size (k, d) 
//
void grad_xi(vector<double> & x, vector<vector<double> > & grad)
{
  int idx ;
  for (int i = 0; i < N; i ++)
  {
    fill(grad[i].begin(), grad[i].end(), 0) ;
    for (int j = 0 ; j < N; j ++)
      grad[i][i * N + j] = x[i * N + j] ;
  }
  idx = N ; 
  for (int i = 0 ; i < N ; i ++)
    for (int j = i+1 ; j < N; j ++)
    {
      fill(grad[idx].begin(), grad[idx].end(), 0) ;
      for (int j0 = 0 ; j0 < N; j0 ++)
      {
	grad[idx][i * N + j0] = x[j * N + j0] ;
	grad[idx][j * N + j0] = x[i * N + j0] ;
      }
      idx ++;
    }
}

void allocate_mem()
{
  qr_mat = (double *) malloc( (d * k) * sizeof(double) ) ;
  qr_work = (double *) malloc( k * sizeof(double) ) ;
  qr_tau = (double *) malloc( k * sizeof(double) ) ;
  mat_a = (double *) malloc( (k * k) * sizeof(double) ) ;
  ipiv = (int *) malloc( k * sizeof(int) ) ;
  linear_sol = (double *) malloc( k * sizeof(double) ) ;
}

void deallocate_meme() 
{
  free(qr_mat) ;
  free(qr_work) ;
  free(qr_tau) ;
  free(mat_a) ;
  free(ipiv) ;
  free(linear_sol) ;
}

void qr_decomp(vector<vector<double> > & grad, vector<vector<double> > & tangent_vec)
{
  int lda , lwork, info ;
  lda = d ;

  for (int i = 0 ; i < d ; i ++)
    for (int j = 0 ; j < k; j ++)
      qr_mat[i * k + j] = grad[j][i] ;
      
  dgeqrf_(&d, &k, qr_mat, &lda, qr_tau, qr_work, &lwork, &info) ;

  if (info != 0)
  {
    printf("return value of QR (step 1) is wrong: info=%d!\n", info) ;
    exit(1) ;
  }

  dorgqr_(&d, &k, &k, qr_mat, &d, qr_tau, qr_work, &lwork, &info) ;

  if (info != 0)
  {
    printf("return value of QR (step 2) is wrong: info=%d\n", info) ;
    exit(1) ;
  }
}

void linear_solver(vector<vector<double> > & mat, vector<double> & sol)
{
  int nrhs, info , lda, ldb ;
  nrhs = 1; lda = k ; ldb = k ; 
  dgesv_(&k, &nrhs, mat_a, &lda, ipiv, linear_sol, &ldb, &info) ;

  extern void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *) ;

  if (info != 0)
  {
    printf("return value of DGESV is wrong: info=%d\n", info) ;
    exit(1) ;
  }

}

double vec_dot(vector<double> & v1, vector<double> & v2) 
{
  double s ;
  s = 0 ;
  for (int i = 0 ; i < d; i ++)
    s += v1[i] * v2[i] ;
  return s;
}

/*
 * Compute the projection by Newton-method.
 *
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
 */

/*
 * update state by MCMC
 *
void update_state(vector<double> & state )
{
  double r, norm , r1, accept_prob , tmp ;
  int flag ;
  vector<double> y_state, yy_state ;

  y_state.resize(2) ;
  yy_state.resize(2) ;
  
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

   * Step 2: projection by Newton-method.
  flag = newton_projection(tmp_state, grad_vec, y_state) ;

  if (flag == 0) {
    forward_newton_counter ++ ;
    return ;
  }

   * normal direction at y 
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
*/

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
  h = 1.00 ;
  T = n * h ;

  // compute the total steps
  output_every_step = 1000 ;

  beta = 1.0 ;

//  eps_tol = 1e-7 ;
  eps_tol = 1e-12 ;
  newton_grad_tol = 1e-7 ; 
  newton_max_step = 10 ;
  reverse_tol = 1e-10 ;

  tot_newton_step = 0 ; 
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

  printf("SO(N), N=");
  cin >> N ;
  d = N * N ;
  k = N * (N+1) / 2 ;

  //initial state to be the identity matrix
  state.resize(d, 0) ;
  for (int i = 0 ; i < N; i++)
    state[i * N] = 1.0 ;

  start = clock() ;

  printf("Total time = %.2f\nh=%.2e\nn=%.1e\nNo. of output states=%d\n", T, h, n *1.0, n / output_every_step) ;
  printf("SO(%d),\t\td=%d\t\tk=%d\n", N, d, k) ;

  for (int i = 0 ; i < n ; i ++)
  {
    //update_state(state) ;
  }

  int tot_rej ;
  out_file.close() ;
  printf("\naverage newton iteration steps = %.2f\n", tot_newton_step * 1.0 / n) ;
  printf("\naverage xi distance = %.3e\n", mean_xi_distance * 1.0 / n) ;
  tot_rej = forward_newton_counter + backward_newton_counter + reverse_check_counter + metropolis_counter ;
  printf("\nRejection rate: Forward\tReverse\tReversibility\tMetrolis\tTotal\n") ;
  printf("\t\t %.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", forward_newton_counter * 1.0 / n, backward_newton_counter * 1.0 / n, reverse_check_counter *1.0/n, metropolis_counter * 1.0 / n, tot_rej * 1.0 / n) ;

  end = clock() ;
  printf("\n\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;

  return 0; 
}
