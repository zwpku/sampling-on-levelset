#include "ex4.h"

int dt_too_large_counter, tot_ode_step ;
double beta, noise_coeff ;

int read_config() ;

vector<double> state, new_state ;

/*
 * vector field of the flow map
 */
void f(vector<double> & x, vector<double> & force)
{
  double s ;
  vector<double> tmp_xi ;
  tmp_xi.resize(k) ;
  xi(x, tmp_xi) ;
  grad_xi(x, grad_vec) ;

  for (int i = 0 ; i < d; i ++)
  {
    force[i] = 0 ;
    for (int j = 0; j < k ; j ++)
      force[i] += tmp_xi[j] * grad_vec[j][i]; 
    force[i] *= -1 ;
  }
}

/*
 * Compute the projection map $\Theta(state)$.
 *
 * The ODE is solved using the Runge-Kutta method. 
 *
 */
int theta_projection(vector<double> & state)
{
  double eps_old, eps, cdt ;
  int step, decrease_counter ;
  vector<double> x0, x, k1, k2, k3 ;

  k1.resize(d) ; k2.resize(d) ; k3.resize(d) ;

  eps = distance_to_levelset(state) ;
  step = 0 ; 
  cdt = dt ;
  decrease_counter = 0 ;
  while (eps > eps_tol) // Runge-Kutta method
  {
    eps_old = eps ;
    x0 = state ; x = state ;
    f(x, k1) ;
    for (int i = 0 ; i < d; i ++)
      x[i] = state[i] + 0.5 * cdt * k1[i] ;
    f(x, k2) ;
    for (int i = 0 ; i < d; i ++)
      x[i] = state[i] + 0.75 * cdt * k2[i] ;
    f(x, k3) ;
    for (int i = 0 ; i < d; i ++)
      state[i] += cdt * (2.0 / 9 * k1[i] + 1.0 / 3 * k2[i] + 4.0 / 9 * k3[i]) ;

    step ++ ;
    eps = distance_to_levelset(state) ;
    if (eps > eps_old)
    {
      state = x0;
      if (verbose_flag == 1)
        log_file << "Warning: distance increased from " << eps_old << " to " << eps << ". Decrease step-size: " << cdt << " to " << cdt * 0.5 << endl ;
      cdt *= 0.5 ;
      eps = eps_old ;
      decrease_counter ++ ;
    }
  }

  tot_ode_step += step ;

  if (decrease_counter > 0)
  {
   if (verbose_flag == 1)
     log_file << "decrease_counter = " << decrease_counter << " (dt, cdt) =" << dt << ' ' << cdt << endl ; 

    dt_too_large_counter ++ ;
  }

  return step ;
}

/*
 * Numerical scheme using the projection map $\Theta$
 *
 */
void update_state_flow_map(vector<double> & state, vector<double> & new_state )
{
  double tmp, r ;
  int step ;

  /* 
   * Step 1: update the state 
   * In our example, matrix a = id.
   */

  for (int i = 0 ; i < d ; i ++)
  {
    r = gennor(0, 1) ;
    new_state[i] = state[i] + noise_coeff * r ;
  }

  tmp = distance_to_levelset(new_state) ;
  mean_xi_distance += tmp ;

  if (verbose_flag == 1)
    log_file << "1. Update : distance = " << tmp << endl << "2. Projection starts:\n" ;

  /* 
   * Step 2: projection using the map $\Theta$.
   */
  step = theta_projection(new_state) ;

  if (verbose_flag == 1)
    log_file << "\tODE steps = " << step << endl << "3. Distance after projection = " << distance_to_levelset(new_state) << endl ;
}

void allocate_mem()
{
  ipiv = (int *) malloc( k * sizeof(int) ) ;
  mat_x = (double *) malloc( d * sizeof(double) ) ;

  /* 
   * gradient vector of xi
   * dimension: k * d
   */
  grad_vec.resize(k) ;
  for (int i = 0 ; i < k; i ++)
    grad_vec[i].resize(d) ;

  trace_series.resize(n) ;
}

void deallocate_mem() 
{
  free(ipiv) ;
  free(mat_x) ;
}

int main ( int argc, char * argv[] ) 
{
  char buf[50] ;
  ofstream out_file ;
  int idx, determinant_flag, new_state_sucess_flag, try_number, tot_try_number ;

  clock_t start , end ;

  read_config() ;

  // step-size in solving ODE or optimization 
  //
  mcmc_flag = 0 ;
  beta = 1.0 ;
  mean_xi_distance = 0 ;
  determinant_tol = 1e-3 ;

  tot_ode_step = 0 ; 
  tot_try_number = 0 ;
  dt_too_large_counter = 0 ;

  init_rand_generator();
  noise_coeff = sqrt(2.0 / beta * h) ;

  printf("SO(N), N=") ;
  cin >> N ;
  d = N * N ;
  k = N * (N+1) / 2 ;

  //start from the identity matrix
  state.resize(d, 0) ;
  for (int i = 0 ; i < N; i++)
    state[i * N + i] = 1.0 ;

  new_state = state ;

  // get memory for several global variables
  allocate_mem() ;

  sprintf(buf, "./data/log_ex4_no_mcmc_son_%d.txt", N) ;
  log_file.open(buf) ;

  start = clock() ;

  printf("\nSO(%d),\td=%d\tk=%d\n", N, d, k) ;
  printf("n=%d\t h=%.3e\toutput_step =%d \n", n, h, output_every_step ) ;

  int progress ;
  progress = 1 ;

  for (int i = 0 ; i < n ; i ++)
  {
    //compute histgram during the simulation
    trace_series[i] = trace(state) ;

    if (i >= progress * 0.01 * n)
    {
      end = clock() ;
      printf("\ni=%d, n=%d, %.1f\%% finished. \n", i, n, progress * 1.0 ) ;
      printf("  %.1fsec, %.1esec per step,  remaining time: %4.1fmin\n  average ODE step=%.2f,  average xi=%.2e\n", (end-start) * 1.0 / CLOCKS_PER_SEC, (end - start) * 1.0 / CLOCKS_PER_SEC / i, (end - start) * 1.0 / CLOCKS_PER_SEC * (n-i) / (i * 60.0), tot_ode_step * 1.0 / i, mean_xi_distance / i) ;
      progress ++ ;
    }

    try_number = 0 ;
    new_state_sucess_flag = 0;
    while (try_number < maximal_try_number)
    {
      update_state_flow_map(state, new_state) ;

      // check the determinant of the new state
      determinant_flag = check_determinant(new_state) ;

      try_number ++ ;

      if (determinant_flag == 0)
      {
	if (verbose_flag == 1)
	  log_file << try_number << " Try: the determinant of new state is not 1! " << endl ;
      } else 
      {
	new_state_sucess_flag = 1;
	break ;
      }
    }

    if (new_state_sucess_flag == 1)
    {
      state = new_state ;
      tot_try_number += try_number ;
      if (verbose_flag == 1)
	  log_file << "\tfound new state after " << try_number << " trials\n" ;
    }
    else {
	if (verbose_flag == 1)
	  log_file << " Error: exceed " << maximal_try_number << " times, still didn't find new state!\n\tMay be the step-size is too large. Exit!" << endl ;
        printf(" Error: exceed %d times, still didn't find new state!\n\t May be the step-size (h=%.4f) is too large. Exit!\n", try_number, h) ;
	exit(1) ;
    }
  }

  printf("\naverage iteration steps = %.2f\n", tot_ode_step * 1.0 / n) ;
  printf("\naverage xi distance = %.3e\n", mean_xi_distance * 1.0 / n) ;
  printf("\nAverage try number = %.1f\n", tot_try_number * 1.0 / n) ;
  printf("\nProb. when dt is large = %.4f\n", dt_too_large_counter * 1.0 / n) ;

  analysis_data_and_output() ;

  log_file.close() ;

  end = clock() ;
  printf("\n\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;

  deallocate_mem() ;

  return 0; 
}
