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

double h, beta , noise_coeff, dt, eps_tol ;
double c, c2, c4 ;
double bin_width, AA ;

vector<double> state, counter_of_each_bin ;

double xi(vector<double> & x)
{
  return 0.5 * (x[0] * x[0] / c2 + x[1] * x[1] - 1.0) ;
}

// non reversible matrix:
// A =  (0 , AA)
//      (-AA, 0)
// 
void f(vector<double> & x, vector<double> & force)
{
  double tmp ;
  tmp = xi(x) ;
  force[0] = -tmp * (x[0] / c2 - AA * x[1]) ;
  force[1] = -tmp * (AA * x[0] / c2 + x[1]) ;
}

void theta_projection(vector<double> & state)
{
  double eps, cdt ;
  int step ;
  vector<double> x, k1, k2, k3 ;

  k1.resize(2) ; k2.resize(2) ; k3.resize(2) ;

  eps = fabs(xi(state)) ;
  step = 0 ; 
  cdt = dt ;
  while (eps > eps_tol) // Runge-Kutta method
  {
    x = state ;
    f(x, k1) ;
    x[0] = state[0] + 0.5 * cdt * k1[0] ; x[1] = state[1] + 0.5 * cdt * k1[1] ;
    f(x, k2) ;
    x[0] = state[0] + 0.75 * cdt * k2[0] ; x[1] = state[1] + 0.75 * cdt * k2[1] ;
    f(x, k3) ;
    state[0] += cdt * (2.0 / 9 * k1[0] + 1.0 / 3 * k2[0] + 4.0 / 9 * k3[0]) ;
    state[1] += cdt * (2.0 / 9 * k1[1] + 1.0 / 3 * k2[1] + 4.0 / 9 * k3[1]) ;
    step ++ ;
    eps = fabs(xi(state)) ;
    // slightly increase the time step-size 
    cdt *= 1.05 ;
  }
  tot_step += step ;
}

void update_state_flow_map(vector<double> & state )
{
  double r ;
  for (int i = 0 ; i < 2 ; i ++)
  {
    r = gennor(0, 1) ;
    state[i] = state[i] + noise_coeff * r ;
  }
  theta_projection(state) ;
}

void update_state_orthogonal_projection(vector<double> & x, double prev_angle )
{
  double angle, df, eps ;
  int step ;
  vector<double> xx ;

  double r ;
  for (int i = 0 ; i < 2 ; i ++)
  {
    r = gennor(0, 1) ;
    x[i] = x[i] + noise_coeff * r ;
  }

  angle = prev_angle ;
  df = (1 - c2) * 0.5 * sin(2 * angle) + c * x[0] * sin(angle) - x[1] * cos(angle) ;

  eps = fabs(df) ;
  step = 0 ; 
  xx.resize(2) ;
  while (eps > eps_tol)
  {
    angle -= df * dt ;
    df = (1 - c2) * 0.5 * sin(2 * angle) + c * x[0] * sin(angle) - x[1] * cos(angle) ;
    xx[0] = c * cos(angle) ; xx[1] = sin(angle) ;

    // check convergence based on 1) derivative; 2) orthogonal condition
    eps = max( fabs(df), fabs(xx[1] * (x[0]-xx[0]) - xx[0]/c2 * (x[1] - xx[1])) ) ;
    step ++ ;
  }

  x = xx; tot_step += step ;
}

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

  switch (scheme_id) {
    case 0 :
      T = 50000 ;
      h = 0.01 ;
      break ;
   case 1 :
      T = 30000 ;
      h = 0.002 ;
      break ;
    case 2 :
      T = 50000 ;
      h = 0.01 ;
      break ;
    case 3 :
      T = 1000 ;
      h = 0.0001 ;
      break ;
  }
  n = int(T / h) ;
  output_every_step = 1000 ;
  dt = 0.1 ;
  beta = 1.0 ;
  c = 3.0 ;
  c2 = c*c ;
  c4 = c2 * c2 ;
  eps_tol = 1e-7 ;
  tot_step = 0 ; 
  n_bins = 50 ;
  bin_width = 2 * pi / n_bins ;
  counter_of_each_bin.resize(n_bins,0) ;

  init_rand_generator();
  noise_coeff = sqrt(2.0 / beta * h) ;

  //initial state
  state.resize(2) ;
  state[0] = c ; state[1] = 0.0 ;

  start = clock() ;

  printf("Scheme=%d\nTotal time = %.2f\nh=%.2e\nn=%.1e\nNo. of output states=%d\n", scheme_id, T, h, n *1.0, n / output_every_step) ;

  sprintf(buf, "../data/traj_%d.txt", scheme_id) ;

  out_file.open(buf) ;

  out_file << n / output_every_step << endl ;

  angle = 0 ;
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
    angle = atan2(state[1], state[0] / c) ;
    // change angle to [0, 2*pi]
    if (angle < 0) angle = 2 * pi + angle ;
    idx = int (angle / bin_width) ;
    counter_of_each_bin[idx] ++ ;
  }

  out_file.close() ;

  printf("\naverage iteration steps = %.2f\n", tot_step * 1.0 / n) ;

  sprintf(buf, "../data/counter_%d.txt", scheme_id) ;
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
