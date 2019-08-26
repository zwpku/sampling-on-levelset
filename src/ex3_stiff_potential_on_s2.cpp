#include "ex3.h"

const double pi = atan(1) * 4 ;

int n, id_mat_a_flag, large_step_size_flag ;
int tot_step ;
int n_bins_theta, n_bins_phi ;

double h, beta , noise_coeff, stiff_eps ;

double bin_width_theta, bin_width_phi ;
double mean_xi_distance ;

vector<double> state, theta_counter_of_each_bin, phi_counter_of_each_bin;

int read_config() ;

// the reaction coordinate function 
// In this example, the level set is the unit sphere S^2 in R^3
double xi(vector<double> & x)
{
  return 0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1.0) ;
}

/*
 * vector field of the flow map, when the matrix a is non-identity
 */
void f(vector<double> & x, vector<double> & force)
{
  double tmp ;
  tmp = xi(x) ;
  force[0] = -tmp * (tmp * 2 + 1) * x[0] ;
  force[1] = -tmp * (tmp * 2 + 1) * x[1] ;
  force[2] = -tmp * (tmp * 2 + 1) * x[2] ;
}

double compute_theta_angle(vector<double> & x)
{
  return asin(x[2]) ;
}

// compute the gradient of the angle theta
void compute_grad_theta(vector<double> & x, vector<double> & grad)
{
  double tmp ;
  tmp = sqrt(x[0]*x[0] + x[1]*x[1]) ;
  grad[0] = -x[0] * x[2] / tmp ;
  grad[1] = -x[1] * x[2] / tmp ;
  grad[2] = tmp ;
}

// return angle \phi \in [0, 2\pi] 
double compute_phi_angle(vector<double> & x)
{
  double tmp ;
  tmp = atan2(x[1], x[0]) ;
  if (tmp < 0) tmp = tmp + 2 * pi ;
  return tmp ;
}

double l2_norm(vector<double> & x)
{
  return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) ;
}

/*
 * Compute the projection map $\Theta(state)$.
 *
 * The ODE is solved using the Runge-Kutta method. 
 *
 */
void theta_projection(vector<double> & state)
{
  double eps, cdt, eps_old ;
  int step ;
  vector<double> x, k1, k2, k3 ;

 // In fact, in both case, the projections by flow maps are simply the orthogonal projection on the sphere.
  /*
  if (id_mat_a_flag == 0) 
  {
    k1.resize(3) ; k2.resize(3) ; k3.resize(3) ;
    eps = fabs(xi(state)) ;
    step = 0 ; 
    cdt = dt ;
    while (eps > eps_tol) // Runge-Kutta method
    {
      eps_old = eps ;
      x = state ;
      f(x, k1) ;
      x[0] = state[0] + 0.5 * cdt * k1[0] ; 
      x[1] = state[1] + 0.5 * cdt * k1[1] ;
      x[2] = state[2] + 0.5 * cdt * k1[2] ;
      f(x, k2) ;
      x[0] = state[0] + 0.75 * cdt * k2[0] ; 
      x[1] = state[1] + 0.75 * cdt * k2[1] ;
      x[2] = state[2] + 0.75 * cdt * k2[2] ;
      f(x, k3) ;
      x = state ;
      state[0] += cdt * (2.0 / 9 * k1[0] + 1.0 / 3 * k2[0] + 4.0 / 9 * k3[0]) ;
      state[1] += cdt * (2.0 / 9 * k1[1] + 1.0 / 3 * k2[1] + 4.0 / 9 * k3[1]) ;
      state[2] += cdt * (2.0 / 9 * k1[2] + 1.0 / 3 * k2[2] + 4.0 / 9 * k3[2]) ;
      step ++ ;
      eps = fabs(xi(state)) ;
      if (eps > eps_old) // if the state becomes further, we decrease the step-size
      {
	printf("eps_old = %.3e\t eps=%.3e\n", eps_old, eps) ;
	state = x ; 
	cdt *= 0.95;
	eps = eps_old ;
      }
//      printf("\tode step=%d, \teps=%.4e, \tstate: (%.4f, %.4f, %.4f)\n", step, xi(state), state[0], state[1], state[2]) ;
    }
    tot_step += step ;
  } 
  */

    double norm ;
    x = state ;
    norm = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) ;
    for (int i = 0 ; i < 3 ; i++)
      state[i] = x[i] / norm  ;
}

void compute_matrix_dadx_sigma(vector<double> & x, vector<vector<double> > & sigma, vector<double> & da)
{
  double norm ;
  sigma.resize(3) ; da.resize(3) ;
  for (int i = 0 ; i < 3 ; i ++)
    sigma[i].resize(3) ;

  norm = sqrt(x[0] * x[0] + x[1] * x[1]) ;

  sigma[0][0] = x[0] ; sigma[0][1] = x[1] ; sigma[0][2] = - sqrt(stiff_eps) * x[0] * x[2] / norm ;
  sigma[1][0] = x[1] ; sigma[1][1] = -x[0] ; sigma[1][2] = - sqrt(stiff_eps) * x[1] * x[2] / norm ;
  sigma[2][0] = x[2] ; sigma[2][1] = 0.0 ; sigma[2][2] = sqrt(stiff_eps) * norm ;

  da[0] = stiff_eps * x[0] * x[2] * x[2] / (norm * norm) + (3 - stiff_eps) * x[0] ;
  da[1] = stiff_eps * x[1] * x[2] * x[2] / (norm * norm) + (3 - stiff_eps) * x[1] ;
  da[2] = (4 - 2 * stiff_eps) * x[2] ;
}

/*
 * Numerical scheme using the projection map $\Theta$
 *
 */
void update_state_flow_map(vector<double> & state )
{
  double theta, tmp ;
  vector<double> grad ;

  grad.resize(3) ;
  compute_grad_theta(state, grad) ;
  theta = compute_theta_angle(state) ;

  /* 
   * Step 1: update the state 
   */
  if (id_mat_a_flag == 1) // if a is an identity matrix
  {
    double r ;
    for (int i = 0 ; i < 3 ; i ++)
    {
      r = gennor(0, 1) ;
      state[i] = state[i] - 1.0 / stiff_eps * h * theta * grad[i] + noise_coeff * r ;
    }
//    printf("norm of grad=%.4e\t\n", 1.0 / stiff_eps * theta * l2_norm(grad) );
  } else 
  {
    vector<double> r, dadx ;
    vector<vector<double> > sigma ;
    r.resize(3) ;

    compute_matrix_dadx_sigma(state, sigma, dadx) ;

    for (int i = 0 ; i < 3 ; i ++)
      r[i] = gennor(0,1) ;
    for (int i = 0 ; i < 3 ; i ++)
    {
      tmp = 0;
      for (int j = 0; j < 3 ; j ++)
	tmp += sigma[i][j] * r[j] ;
      state[i] = state[i] + (-theta * grad[i] + 1.0 / beta * dadx[i]) * h + noise_coeff * tmp ;
    }
//    printf("norm of grad=%.4e\t norm of dadx=%.4e\t\n", theta * l2_norm(grad), l2_norm(dadx) );

//    printf("\tafter update: state=(%.4f, %.4f, %.4f)\n", state[0], state[1], state[2]) ;
  }

  mean_xi_distance += fabs(xi(state)) ;
  /* 
   * Step 2: projection using the map $\Theta$.
   */
  theta_projection(state) ;
//  printf("\tafter projection state=(%.4f, %.4f, %.4f)\n", state[0], state[1], state[2]) ;
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
  double theta_angle, phi_angle, T;

  clock_t start , end ;

  // initialize

  read_config() ;

  beta = 1.0 ;

  printf("stiff_eps=%.4e\n", stiff_eps) ;

  if (id_mat_a_flag == 0)
  {
    h = 0.01 ;
    printf("non-identity matrix, h = %.4f\n", h) ;
  }
  else {
    if (large_step_size_flag == 1)
      h = 0.005 ;
    else h=0.0002 ;

    printf("a=id, h = %.4f\n", h) ;
  }

  T = n * h ;
  output_every_step = 1000 ;
  mean_xi_distance = 0 ;

  tot_step = 0 ; 
  n_bins_theta = 500 ;
  n_bins_phi = 100 ;

  // divied [-pi/2, pi/2] to n_bins_theta with equal width
  bin_width_theta = pi / n_bins_theta ;
  // divied [0, 2pi] to n_bins_phi with equal width
  bin_width_phi = 2 * pi / n_bins_phi ;

  theta_counter_of_each_bin.resize(n_bins_theta,0) ;
  phi_counter_of_each_bin.resize(n_bins_phi,0) ;

  init_rand_generator();
  noise_coeff = sqrt(2.0 / beta * h) ;

  //initial state
  state.resize(3) ;
  state[0] = 1.0 ; state[1] = 0.0 ; state[2] = 0.0 ;

  start = clock() ;

  printf("n=%d\t\th=%.2e\tT=%.2f\t\tNo. of output states=%d\n", n, h, T, n / output_every_step) ;

  sprintf(buf, "./data/ex3_traj_%d_%d.txt", id_mat_a_flag, large_step_size_flag) ;

  out_file.open(buf) ;

  out_file << n / output_every_step << endl ;

  for (int i = 0 ; i < n ; i ++)
  {
    if (i % output_every_step == 0)
    {
      out_file << state[0] << ' ' << state[1] << ' ' << state[2] << endl ;
    }

    update_state_flow_map(state) ;

    //compute histgram of angle theta and phi during the simulation

    // phi in [0, 2pi]
    phi_angle = compute_phi_angle(state) ;
    idx = int (phi_angle / bin_width_phi) ;
    phi_counter_of_each_bin[idx] ++ ;

    // theta in [-pi/2, pi/2]
    theta_angle = compute_theta_angle(state) ;
    idx = int ((theta_angle + pi * 0.5) / bin_width_theta) ;
    theta_counter_of_each_bin[idx] ++ ;
  }

  out_file.close() ;

  printf("\naverage iteration steps = %.2f\n", tot_step * 1.0 / n) ;
  printf("\naverage xi distance = %.3e\n", mean_xi_distance * 1.0 / n) ;

  sprintf(buf, "./data/ex3_theta_counter_%d_%d.txt", id_mat_a_flag, large_step_size_flag) ;
  out_file.open(buf) ;
  out_file << n << ' ' << n_bins_theta << endl ;
  for (int i = 0 ; i < n_bins_theta ; i++)
  {
    out_file << theta_counter_of_each_bin[i] << ' ';
  }
  out_file << endl ;
  out_file.close() ;

  sprintf(buf, "./data/ex3_phi_counter_%d_%d.txt", id_mat_a_flag, large_step_size_flag) ;
  out_file.open(buf) ;
  out_file << n << ' ' << n_bins_phi << endl ;
  for (int i = 0 ; i < n_bins_phi ; i++)
  {
    out_file << phi_counter_of_each_bin[i] << ' ';
  }
  out_file << endl ;
  out_file.close() ;

  end = clock() ;
  printf("\n\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;

  return 0; 
}
