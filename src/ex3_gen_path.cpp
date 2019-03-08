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

int n, scheme_id ;
int tot_step ;
int n_bins ;

int ref_ode_sol_flag, ref_ode_tot_step, ref_ode_c_step ;
double ref_ode_distance, ode_rescale_factor ;

double h, beta , noise_coeff, dt, eps_tol ;

/*
 * c2=c^2, c4=c^4
 */
double c, c2, c4 ;
double bin_width, AA ;

double mean_xi_distance ;
double pot_coeff ; 

vector<double> state, counter_of_each_bin ;

ofstream log_file ;

// the reaction coordinate function 
double xi(vector<double> & x)
{
  return 0.5 * (x[0] * x[0] / c2 + x[1] * x[1] - 1.0) ;
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

// non reversible matrix:
// A =  (0 , AA)
//      (-AA, 0)
// 

/*
 * vector field of the flow map
 */
void f(vector<double> & x, vector<double> & force, int ref_flag )
{
  double tmp ;
  tmp = xi(x) ;
  force[0] = -tmp * (x[0] / c2 - AA * x[1]) ;
  force[1] = -tmp * (AA * x[0] / c2 + x[1]) ;

  if (ref_flag == 0)
  {
    force[0] *=  pow(fabs(tmp), -ode_rescale_factor) ;
    force[1] *=  pow(fabs(tmp), -ode_rescale_factor) ;
  }
}

double norm(vector<double> & x)
{
  return sqrt(x[0] * x[0] + x[1] * x[1]) ;
}

/*
 * Compute the projection map $\Theta(state)$.
 *
 * The ODE is solved using the Runge-Kutta method. 
 *
 */
int theta_projection( vector<double> & state, double init_dt, int ref_flag )
{
  double eps, cdt , eps_old ;
  int step ;
  vector<double> x, x0, k1, k2, k3 ;

  k1.resize(2) ; k2.resize(2) ; k3.resize(2) ;

  eps = fabs(xi(state)) ;
  step = 0 ; 
  cdt = init_dt ;
  while (eps > eps_tol) // Runge-Kutta method
  {
    eps_old = eps ;
    x = state ; x0 = state ;
    f(x, k1, ref_flag) ;
    x[0] = state[0] + 0.5 * cdt * k1[0] ; x[1] = state[1] + 0.5 * cdt * k1[1] ;
    f(x, k2, ref_flag) ;
    x[0] = state[0] + 0.75 * cdt * k2[0] ; x[1] = state[1] + 0.75 * cdt * k2[1] ;
    f(x, k3, ref_flag) ;
    state[0] += cdt * (2.0 / 9 * k1[0] + 1.0 / 3 * k2[0] + 4.0 / 9 * k3[0]) ;
    state[1] += cdt * (2.0 / 9 * k1[1] + 1.0 / 3 * k2[1] + 4.0 / 9 * k3[1]) ;
    step ++ ;
    eps = fabs(xi(state)) ;

    if (eps > eps_old)
    {
      state = x0 ;
      //log_file << "Warning: in step " << step << ", distance increased from " << eps_old << " to " << eps << ". Decrease step-size: " << cdt << " to " << cdt * 0.5 << endl ;
      cdt *= 0.5 ;
      eps = eps_old ;
    }
  }
  if (ref_flag == 0) tot_step += step ;

  return step ;
}

/*
 * Numerical scheme using the projection map $\Theta$
 *
 */
void update_state_flow_map(vector<double> & state )
{
  vector<double> grad ;
  grad.resize(2) ;

  double r ;
  /* 
   * Step 1: update the state 
   * Notice that, in our example, matrix a = id.
   */

  grad_U(state, grad) ;
  for (int i = 0 ; i < 2 ; i ++)
  {
    r = gennor(0, 1) ;
    state[i] = state[i] - grad[i] * h + noise_coeff * r ;
  }

  mean_xi_distance += fabs(xi(state)) ;

  // test accuracy using a small step-size

  vector<double> new_state_tmp ;
  if (ref_ode_sol_flag == 1) new_state_tmp = state ;

  /* 
   * Step 2: projection using the map $\Theta$.
   */
  int step ;
  step = theta_projection(state, dt, 0) ;

  if (ref_ode_sol_flag == 1)  // compare to the reference solution of the ODE
  {
    int step_test ;
    // solve the ODE with small step-size, no rescaling.
    step_test = theta_projection(new_state_tmp, dt * 0.01, 1) ;

    // compute the l^2 distance of the ODE solutions.
    double tmp ;
    tmp = 0;
    for (int ii = 0; ii < 2; ii++)
      tmp += fabs(state[ii] - new_state_tmp[ii]) * fabs(state[ii] - new_state_tmp[ii]) ;

    ref_ode_distance += sqrt(tmp) ;

    log_file << "Distance of ODE solutions=" << sqrt(tmp) << ", ODE steps =" << step << ',' << step_test << endl ; 

    ref_ode_c_step ++ ;
  }
}

/*
 * Scheme using projection along geodesic curves
 */

void update_state_orthogonal_projection(vector<double> & x, double prev_angle )
{
  double angle, df, eps ;
  int step ;
  vector<double> xx ;

  double r ;

  /*
   * Step 1: update the state
   */
  for (int i = 0 ; i < 2 ; i ++)
  {
    r = gennor(0, 1) ;
    x[i] = x[i] + noise_coeff * r ;
  }

  /*
   * Solve the optimization problem, parametrized by the angle theta
   */

  // initialize using the solution from the previous step
  angle = prev_angle ;

  // derivative of the objective function 
  df = (1 - c2) * 0.5 * sin(2 * angle) + c * x[0] * sin(angle) - x[1] * cos(angle) ;

  eps = fabs(df) ;
  step = 0 ; 
  xx.resize(2) ;
  while (eps > eps_tol)
  {
    // gradient descent
    angle -= df * dt ;
    // compute the derivative of the objective function 
    df = (1 - c2) * 0.5 * sin(2 * angle) + c * x[0] * sin(angle) - x[1] * cos(angle) ;
    xx[0] = c * cos(angle) ; xx[1] = sin(angle) ;

    /* 
     * Check convergence based on 
     * 1) derivative; 2) orthogonal condition
     */
    eps = max( fabs(df), fabs(xx[1] * (x[0]-xx[0]) - xx[0]/c2 * (x[1] - xx[1])) ) ;
    step ++ ;
  }

  x = xx; tot_step += step ;
}

/*
 * Euler-Maruyama discretization of the SDE
 */
void update_state_euler_sde(vector<double> & x )
{
  double r1, r2 ;
  double tmp1, tmp2, p11, p12, p22 ;
  r1 = gennor(0, 1) ; r2 = gennor(0, 1) ;
  p11 = x[1] * x[1] / (x[0] * x[0] / c4 + x[1] * x[1]) ;
  p12 = - c2 * x[0] * x[1] / (c4 * x[1] * x[1] + x[0] * x[0]) ;
  p22 = x[0] * x[0] / (x[0] * x[0] + c4 * x[1] * x[1]) ;
  tmp1 = x[0] * ((c2-2) * x[1] * x[1] - x[0] * x[0] / c2) / pow(c2 * x[1] * x[1] + x[0] * x[0] / c2 , 2);
  tmp2 = x[1] * (-c2 * x[1] * x[1] + x[0] * x[0] * (1.0 / c2 - 2)) / pow(c2 * x[1] * x[1] + x[0] * x[0] / c2, 2);
  x[0] += tmp1 / beta * h + noise_coeff * (p11 * r1 + p12 * r2) ;
  x[1] += tmp2 / beta * h + noise_coeff * (p12 * r1 + p22 * r2) ;
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
  double angle , T;

  clock_t start , end ;

  // initialize

  cout << "0: \\Theta\n1: \\Theta^{A}\n2: \\Pi\n3: E-M\n" << "input scheme id= " ;
  cin >> scheme_id  ;

  /*
   * initialize the total time T and the step-size h.
   */
  switch (scheme_id) {
    case 0 : 
      n = 20000000 ;
      h = 0.010 ;
      break ;
   case 1 :
      n = 20000000;
      h = 0.005 ;
      break ;
    case 2 :
      n = 20000000;
      h = 0.01 ;
      break ;
    case 3 :
      n = 20000000;
      h = 0.0001 ;
      break ;
  }

  // compute the total steps
  T = n * h ;
  output_every_step = 1000 ;
  // step-size in solving ODE or optimization 
  dt = 0.1 ;
  beta = 1.0 ;
  mean_xi_distance = 0 ;
  pot_coeff = 0.0 ;

  // this is the parameter \kappa in the paper
  ode_rescale_factor = 0.5 ;

  c = 3.0 ;
  c2 = c*c ;
  c4 = c2 * c2 ;

//  eps_tol = 1e-7 ;
  eps_tol = 1e-8 ;
  tot_step = 0 ; 
  n_bins = 50 ;

  // divied [0, 2pi] to n_bins with equal width
  bin_width = 2 * pi / n_bins ;

  counter_of_each_bin.resize(n_bins,0) ;

  init_rand_generator();
  noise_coeff = sqrt(2.0 / beta * h) ;

  //initial state
  state.resize(2) ;
  state[0] = c ; state[1] = 0.0 ;

  start = clock() ;

  printf("Scheme=%d\nTotal time = %.2f\nh=%.2e\nn=%.1e\nNo. of output states=%d\n", scheme_id, T, h, n *1.0, n / output_every_step) ;

  sprintf(buf, "./data/ex3_traj_%d.txt", scheme_id) ;

  out_file.open(buf) ;

  out_file << n / output_every_step << endl ;

  sprintf(buf, "./data/log_ex3_scheme_%d.txt", scheme_id) ;
  log_file.open(buf) ;

  angle = 0 ;

  ref_ode_distance = 0.0 ;
  ref_ode_tot_step = 20 ;
  ref_ode_c_step = 0;
  ref_ode_sol_flag = 1;

  for (int i = 0 ; i < n ; i ++)
  {
    if (i % output_every_step == 0)
      out_file << state[0] << ' ' << state[1] << endl ;

    switch (scheme_id) {
      case 0 : 
	AA = 0.0; 
	update_state_flow_map(state) ;
	break ;
      case 1 :
	AA = 0.5; 
	update_state_flow_map(state) ; // the same as in the case scheme_id=0, but with nonzero AA.
	break ;
      case 2 :
	update_state_orthogonal_projection(state, angle) ;
	break ;
      case 3 :
	update_state_euler_sde(state) ;
	break ;
    }

    if ( (ref_ode_c_step == ref_ode_tot_step) && (scheme_id < 2) )
      {
	ref_ode_sol_flag = 0 ;
	ref_ode_c_step ++ ;
	printf("average ODE error = %.2e\n", ref_ode_distance / ref_ode_tot_step) ;
      }


    //compute histgram during the simulation
    angle = atan2(state[1], state[0] / c) ;
    // change angle to [0, 2*pi]
    if (angle < 0) angle = 2 * pi + angle ;
    idx = int (angle / bin_width) ;
    counter_of_each_bin[idx] ++ ;
  }

  out_file.close() ;

  printf("\naverage iteration steps = %.2f\n", tot_step * 1.0 / n) ;
  printf("\naverage xi distance = %.3e\n", mean_xi_distance * 1.0 / n) ;

  sprintf(buf, "./data/ex3_counter_%d.txt", scheme_id) ;
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

  log_file.close() ;

  return 0; 
}
