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

double eps_tol, reverse_tol, newton_grad_tol , size_s ;

double mean_xi_distance ;

int newton_max_step ;

int forward_newton_counter, metropolis_counter, backward_newton_counter, reverse_check_counter ;

double bin_width ;

vector<double> state, v_vec, v_vec_prime, y_state, tmp_state, counter_of_each_bin ;

vector<vector<double> > grad_vec , tangent_vec_array, grad_vec_y, tangent_vec_array_y ;

// arrays for QR decomposition
double * qr_tau, * qr_work , * qr_grad_mat ;

double * mat_a, * linear_sol ; 
int *ipiv ;

// used for debug
double orthogonal_tol ;

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
  qr_grad_mat = (double *) malloc( (d * d) * sizeof(double) ) ;
  qr_work = (double *) malloc( d * sizeof(double) ) ;
  qr_tau = (double *) malloc( k * sizeof(double) ) ;
  mat_a = (double *) malloc( (k * k) * sizeof(double) ) ;
  ipiv = (int *) malloc( k * sizeof(int) ) ;
  linear_sol = (double *) malloc( k * sizeof(double) ) ;

  // (d-k) orthogonal vectors on the tangent space
  tangent_vec_array.resize(d-k) ;
  for (int i = 0 ; i < d - k ; i++)
    tangent_vec_array[i].resize(d) ;

  /* 
   * gradient vector of xi
   * dimension: k * d
   */
  grad_vec.resize(k) ;
  for (int i = 0 ; i < k; i ++)
    grad_vec[i].resize(d) ;

  grad_vec_y = grad_vec ;
  tangent_vec_array_y = tangent_vec_array ;

  v_vec.resize(d) ;
  v_vec_prime.resize(d) ;
}

void deallocate_mem() 
{
  free(qr_grad_mat) ;
  free(qr_work) ;
  free(qr_tau) ;
  free(mat_a) ;
  free(ipiv) ;
  free(linear_sol) ;
}

double vec_dot(vector<double> & v1, vector<double> & v2) 
{
  double s ;
  s = 0 ;
  for (int i = 0 ; i < v1.size(); i ++)
    s += v1[i] * v2[i] ;
  return s;
}

// check orthogonality
int check_orthognality()
{
  double tmp ;
  int flag ;
  printf("checking inner products between tangent vectors ... ") ;
  flag = 1 ;
  for (int i = 0 ; i < d-k; i ++)
    for (int j = i ; j < d-k; j ++)
    {
      tmp = vec_dot(tangent_vec_array[i], tangent_vec_array[j]) ;
      if ( (j == i) && (fabs(tmp - 1) > orthogonal_tol) ) 
      {
	printf("Warning: orthogonality check failed!\n");
	printf("|<v_%d, v_%d> - 1| = %.4e > %.4e\n", i, j, fabs(tmp-1), orthogonal_tol) ;
	flag = 0;
      } 
      if ( (j > i) && (fabs(tmp) > orthogonal_tol) ) 
      {
	printf("\nWarning: orthogonality check failed!\n");
	printf("|<v_%d, v_%d>| = %.4e > %.4e\n", i, j, fabs(tmp), orthogonal_tol) ;
	flag = 0;
      } 
    }

  if (flag == 1) printf("passed\n") ;

  printf("check inner products between tangents and gradient of \\xi ... ") ;
  for (int i = 0 ; i < k; i ++)
    for (int j = 0 ; j < d-k; j ++)
    {
      tmp = vec_dot(grad_vec[i], tangent_vec_array[j]) ;

      if (fabs(tmp) > orthogonal_tol) 
      {
	printf("\nWarning: orthogonality check failed!\n");
	printf("|<n_%d, v_%d>| = %.4e > %.4e\n", i, j, fabs(tmp), orthogonal_tol) ;
	flag = 0;
      } 
    }

  if (flag == 1) printf("passed\n") ;

  return flag ;
}

void qr_decomp(vector<vector<double> > & grad, vector<vector<double> > & t_vec)
{
  int lda , lwork, info ;

  printf("\nQR decomposition...\n") ;

  lda = d ;
  lwork = k ;
  for (int j = 0 ; j < k; j ++)
    for (int i = 0 ; i < d ; i ++)
      qr_grad_mat[j * d + i] = grad[j][i] ;
      
  dgeqrf_(&d, &k, qr_grad_mat, &lda, qr_tau, qr_work, &lwork, &info) ;

  if (info != 0)
  {
    printf("return value of QR (step 1) is wrong: info=%d!\n", info) ;
    exit(1) ;
  }

  lwork = d ;
  dorgqr_(&d, &d, &k, qr_grad_mat, &d, qr_tau, qr_work, &lwork, &info) ;

  if (info != 0)
  {
    printf("return value of QR (step 2) is wrong: info=%d\n", info) ;
    exit(1) ;
  }

  // store the (d-k) tangent vectors
  for (int i = 0 ; i < d-k; i ++)
  {
    for (int j = 0 ; j < d; j ++)
      tangent_vec_array[i][j] = qr_grad_mat[(k+i)*d + j] ;
  }

  printf("QR decomposition finished.\n\n") ;

  check_orthognality() ;
}

void linear_solver(vector<vector<double> > & mat, vector<double> & sol)
{
  int nrhs, info , lda, ldb ;
  nrhs = 1; lda = k ; ldb = k ; 
  dgesv_(&k, &nrhs, mat_a, &lda, ipiv, linear_sol, &ldb, &info) ;

  if (info != 0)
  {
    printf("return value of DGESV is wrong: info=%d\n", info) ;
    exit(1) ;
  }
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

// randomly generate a tangent vector 
void generate_v(vector<double> & v_vec) 
{
  vector<double> r ;
  r.resize(d-k) ;
  for (int i = 0 ; i < d-k; i ++)
    r[i] = gennor(0, size_s) ;

  for (int i = 0 ; i < d; i ++)
  {
    v_vec[i] = 0 ;
    for (int j = 0 ; j < d-k; j ++)
      v_vec[i] += r[j] * tangent_vec_array[j][i] ;
  }
}

int projection_by_Newton(vector<vector<double> > & Qx, vector<double> & x0)
{
  int nrhs, info , lda, ldb, step ;
  double * vec_a, s, eps ;
  vector<double> tmp_x, tmp_rhs ;
  vector<vector<double> > tmp_grad_mat ;

  tmp_rhs.resize(k) ;
  tmp_grad_mat.resize(k) ;
  for (int i =0; i < k ; i ++)
    tmp_grad_mat[i].resize(d); 

  // initial solution is set to zero
  vec_a = (double *) malloc( k * sizeof(double) ) ;
  for (int i = 0; i < k; i ++)
    vec_a[i] = 0 ;

  // initial state
  tmp_x = x0 ;
  xi(tmp_x, tmp_rhs) ;
  eps = sqrt(vec_dot(tmp_rhs, tmp_rhs)) ;
  step = 0 ;

  // Newton iteration method 
  while (eps > eps_tol)
  {
    if (step + 1 > newton_max_step) 
    {
      tot_newton_step += step ;
      return 0 ;
    }

    // calculate the derivative at current state
    grad_xi(tmp_x, tmp_grad_mat) ;

    // initialize the matrix mat_a
    for (int i = 0 ; i < k; i ++)
      for (int j = 0 ; j < k; j ++)
      {
	s = 0 ;
	for (int i0 = 0; i0 < d; i0++)
	  s += tmp_grad_mat[i][i0] * Qx[j][i0] ;
	mat_a[i*k+j] = s ;
      }

    // right hand side
    for (int i = 0 ; i < k; i ++)
      linear_sol[i] = -tmp_rhs[i] ;

    nrhs = 1; lda = k ; ldb = k ; 
    dgesv_(&k, &nrhs, mat_a, &lda, ipiv, linear_sol, &ldb, &info) ;
    if (info != 0)
    {
      printf("Warning: return value of the linear solver DGESV is wrong: info=%d\n", info) ;
      tot_newton_step += step + 1 ;
      return 0 ;
    }

    // update the coefficient vector a
    for (int i =0; i < k; i ++)
      vec_a[i] += linear_sol[i] ;

    step ++ ;
    for (int i =0 ; i < d; i ++)
    {
      s = 0;
      for (int j = 0 ; j < k ; j++)
	s += vec_a[j] * Qx[j][i];
      tmp_x[i] = x0[i] + s;
    }

    xi(tmp_x, tmp_rhs) ;
    eps = sqrt(vec_dot(tmp_rhs, tmp_rhs)) ;
  }

  // Newton iteration is successful if we are here
  x0 = tmp_x ;
  tot_newton_step += step ;

  return 1;
}

void compute_v_prime(vector<double> & state, vector<double> & tmp_state, vector<vector<double> > & tangent_vec_array_y) 
{
  double s ;
  vector<double> state_diff , coeff ;
  state_diff.resize(d) ;
  coeff.resize(d-k) ;

  // x-y
  for (int i = 0; i < d; i ++)
    state_diff[i] = state[i] - tmp_state[i] ;

  // coefficients corresponding to the orthonormal basis
  for (int i = 0 ; i < d-k; i++)
    coeff[i] = vec_dot(state_diff, tangent_vec_array_y[i]) ;

  for (int i =0; i < d; i ++)
  {
    s = 0 ;
    for (int j = 0 ; j < d-k; j ++)
      s+= coeff[j] * tangent_vec_array_y[j][i];
    v_vec_prime[i] = s;
  }
}

int main ( int argc, char * argv[] ) 
{
  char buf[50] ;
  ofstream out_file ;
  int idx ;
  int newton_success_flag ;
  int output_every_step ;
  double angle , T , accept_prob , tmp ;

  clock_t start , end ;

  n = 5000000 ;
  h = 1.00 ;
  T = n * h ;

  // compute the total steps
  output_every_step = 1000 ;

  beta = 1.0 ;

//  eps_tol = 1e-7 ;
  eps_tol = 1e-12 ;
  size_s = 1.0 ;
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

  tmp_state.resize(d) ;
  y_state.resize(d) ;
  allocate_mem() ;

  start = clock() ;

  printf("\nTotal time = %.2f\nh=%.2e\nn=%.1e\nNo. of output states=%d\n", T, h, n *1.0, n / output_every_step) ;
  printf("\nSO(%d),\td=%d\tk=%d\n", N, d, k) ;

  // for the initial state, compute the Jaccobi (gradient) matrix of xi at current state 
  grad_xi(state, grad_vec) ;

  // for the initial state, fine the orthogonal vectors of tangent space by QR decomposition 
  qr_decomp(grad_vec, tangent_vec_array) ;

  for (int i = 0 ; i < n ; i ++)
  {
    // randomly generate a vector on the (d-k)-dimensional tangent space
    generate_v(v_vec) ;

    // move along tangent vector v
    for (int i = 0 ; i < d; i ++)
      tmp_state[i] = state[i] + v_vec[i] ;

    // projection by newton method
    newton_success_flag = projection_by_Newton(grad_vec, tmp_state) ;

    if (newton_success_flag == 0) // increase the counter, if we didn't find a new state
    {
      forward_newton_counter ++ ;
      continue;
    }

    /* 
     * we found a NEW state: tmp_state. Now we do Metropolis-Hasting 
     *
     */

    y_state = tmp_state ;

    // compute the Jaccobi (gradient) matrix of xi at the NEW state 
    grad_xi(y_state, grad_vec_y) ;

    // fine the orthogonal vectors of tangent space (NEW state) by QR decomposition 
    qr_decomp(grad_vec_y, tangent_vec_array_y) ;

    // decide vector v'
    compute_v_prime(state, y_state, tangent_vec_array_y) ;

    // density ratio of two vectors v and v'
    accept_prob = exp(-(vec_dot(v_vec_prime, v_vec_prime) - vec_dot(v_vec, v_vec)) * 0.5 / (size_s * size_s) ) ;
    if (accept_prob > 1.0) accept_prob = 1.0 ;

    tmp = ranf() ;
    if (tmp > accept_prob) // rejected 
    {
      metropolis_counter ++ ;
      continue ;
    }

    /*
     * We have passed the Metropolis-Hasting check. 
     * In the following, we do reversibility check.
     */

    // move state y along tangent vector v'
    for (int i = 0 ; i < d; i ++)
      tmp_state[i] = y_state[i] + v_vec_prime[i] ;

    // projection by newton method
    newton_success_flag = projection_by_Newton(grad_vec_y, tmp_state) ;

    if (newton_success_flag == 0) // increase the counter, if we didn't find a new state
    {
      backward_newton_counter ++ ;
      continue;
    }

    // move to the new state
    state = y_state ;
    // update the gradients and the tangent vectors
    grad_vec = grad_vec_y ;
    tangent_vec_array = tangent_vec_array ;
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

  deallocate_mem() ;

  return 0; 
}
