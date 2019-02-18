#include "ex4_no_mcmc.h"

const double pi = atan(1) * 4 ;

int n, N, d, k, tot_step ;
int n_bins, output_every_step, maximal_try_number, tot_try_number ;

double h, beta, noise_coeff, dt, eps_tol ;

double bin_width, trace_b ;

double mean_xi_distance ;

int verbose_flag ;
ofstream log_file ;

int read_config() ;

vector<double> state, new_state, counter_of_each_bin ;

vector<vector<double> > grad_vec;

double * mat_x ; 
int *ipiv ;

double determinant_tol ;

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

double trace(vector<double> & x)
{
  double s;
  s = 0 ;
  for (int i = 0 ; i < N; i ++)
    s += x[i * N + i];
  return s;
}

int check_determinant(vector<double> & x)
{
  int info, flag ;
  double s ;

  for (int i = 0 ; i < N * N ; i ++)
    mat_x[i] = x[i] ;

  dgetrf_(&N, &N, mat_x, &N, ipiv, &info);

  s = 1.0 ;
  for (int i = 0 ; i < N; i ++)
    s *= mat_x[i * N + i] ;

  for (int i = 0 ; i < N; i ++)
    if (ipiv[i] != i+1) s *= -1 ;

  if (fabs(s-1.0) > determinant_tol)
  {
    flag = 0 ;
    if (verbose_flag == 1)
      log_file << "determinant check failed!\t det=" << std::scientific << s << "\n" ;
  }
  else flag = 1;

  return flag ;
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

double distance_to_levelset(vector<double> & x)
{
  vector<double> tmp_vec ;
  double eps ;
  tmp_vec.resize(k) ;
  xi(x, tmp_vec) ;

  eps = 0 ;
  for (int i = 0 ; i < k; i ++)
    eps += tmp_vec[i] * tmp_vec[i] ;

  return sqrt(eps) ;
}

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
  double eps ;
  int step ;
  vector<double> x, k1, k2, k3 ;

  k1.resize(d) ; k2.resize(d) ; k3.resize(d) ;

  eps = distance_to_levelset(state) ;
  step = 0 ; 
  while (eps > eps_tol) // Runge-Kutta method
  {
    x = state ;
    f(x, k1) ;
    for (int i = 0 ; i < d; i ++)
      x[i] = state[i] + 0.5 * dt * k1[i] ;
    f(x, k2) ;
    for (int i = 0 ; i < d; i ++)
      x[i] = state[i] + 0.75 * dt * k2[i] ;
    f(x, k3) ;
    for (int i = 0 ; i < d; i ++)
      state[i] += dt * (2.0 / 9 * k1[i] + 1.0 / 3 * k2[i] + 4.0 / 9 * k3[i]) ;

    step ++ ;
    eps = distance_to_levelset(state) ;
  }
  tot_step += step ;

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
  int idx, determinant_flag, new_state_sucess_flag, try_number ;
  double T, tmp ;

  clock_t start , end ;

  read_config() ;
  // compute the total steps
  T = n * h ;

  // step-size in solving ODE or optimization 
  //
  beta = 1.0 ;
  mean_xi_distance = 0 ;
  determinant_tol = 1e-3 ;

  tot_step = 0 ; 

  // statistics for output 
  // divied [-trace_b, trace_b] to n_bins with equal width
  bin_width = 2.0 * trace_b  / n_bins ;
  counter_of_each_bin.resize(n_bins, 0) ;

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
  printf("n=%d\t output_step =%d \n", n, output_every_step ) ;

  sprintf(buf, "./data/ex4_no_mcmc_traj_%d.txt", N) ;
  out_file.open(buf) ;
  out_file << n / output_every_step << endl ;

  for (int i = 0 ; i < n ; i ++)
  {
    //compute histgram during the simulation
    tmp = trace(state) ;
    idx = int ((tmp + trace_b) / bin_width) ;
    counter_of_each_bin[idx] ++ ;

    if (i % output_every_step == 0)
      out_file << tmp << ' ' ;

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

  out_file.close() ;

  printf("\naverage iteration steps = %.2f\n", tot_step * 1.0 / n) ;
  printf("\naverage xi distance = %.3e\n", mean_xi_distance * 1.0 / n) ;
  printf("\nAverage try number = %.1f\n", tot_try_number * 1.0 / n) ;

  sprintf(buf, "./data/ex4_no_mcmc_counter_%d.txt", N) ;
  out_file.open(buf) ;

  out_file << n << ' ' << trace_b << ' ' << n_bins << endl ;

  for (int i = 0 ; i < n_bins ; i++)
  {
    out_file << counter_of_each_bin[i] << ' ';
  }
  out_file << endl ;
  out_file.close() ;
  log_file.close() ;

  end = clock() ;
  printf("\n\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;

  deallocate_mem() ;

  return 0; 
}
