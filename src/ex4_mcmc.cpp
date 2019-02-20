#include "ex4.h"

// total steps
int tot_newton_step, tot_success_newton_step ;

double reverse_tol, size_s ;

int newton_solver_counter, newton_solver_converge_counter ;

int forward_newton_counter, metropolis_counter, backward_newton_counter, reverse_check_counter, determinant_counter ;

vector<double> state, v_vec, v_vec_prime, y_state, tmp_state ;

double len_v_bin_width, len_v_bound, max_v_norm_converged, min_v_norm_not_converged ;
int n_v_bins ;
vector<double> bin_len_v_newton_converged, bin_len_v_newton_not_converged ;

vector<vector<double> > tangent_vec_array, grad_vec_y, tangent_vec_array_y ;

// arrays for QR decomposition
double * qr_tau, * qr_work , * qr_grad_mat ;
double * mat_a, * linear_sol ; 

// used for debug
double orthogonal_tol ;

int read_config() ;

void allocate_mem()
{
  qr_grad_mat = (double *) malloc( (d * d) * sizeof(double) ) ;
  qr_work = (double *) malloc( d * sizeof(double) ) ;
  qr_tau = (double *) malloc( k * sizeof(double) ) ;
  mat_a = (double *) malloc( (k * k) * sizeof(double) ) ;
  ipiv = (int *) malloc( k * sizeof(int) ) ;
  linear_sol = (double *) malloc( k * sizeof(double) ) ;
  mat_x = (double *) malloc( d * sizeof(double) ) ;

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
  trace_series.resize(n) ;
}

void deallocate_mem() 
{
  free(qr_grad_mat) ;
  free(qr_work) ;
  free(qr_tau) ;
  free(mat_a) ;
  free(ipiv) ;
  free(linear_sol) ;
  free(mat_x) ;
}

// check orthogonality
int check_orthognality(vector<vector<double> > & grad, vector<vector<double> > & t_vec )
{
  double tmp ;
  int flag ;

  if (verbose_flag == 1)
    log_file << "\tcheck inner products between tangent vectors ... " ;

  flag = 1 ;
  for (int i = 0 ; i < d-k; i ++)
    for (int j = i ; j < d-k; j ++)
    {
      tmp = vec_dot(t_vec[i], t_vec[j]) ;
      if ( (j == i) && (fabs(tmp - 1) > orthogonal_tol) ) 
      {
	printf("\nWarning: orthogonality check failed!\n");
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

  if ( (verbose_flag == 1) && (flag == 1) )
    log_file << "passed" << endl ;

  if (verbose_flag == 1)
    log_file << "\tcheck inner products between tangents and gradient of \\xi ... " ;

  for (int i = 0 ; i < k; i ++)
    for (int j = 0 ; j < d-k; j ++)
    {
      tmp = vec_dot(grad[i], t_vec[j]) ;

      if (fabs(tmp) > orthogonal_tol) 
      {
	printf("\nWarning: orthogonality check failed!\n");
	printf("|<n_%d, v_%d>| = %.4e > %.4e\n", i, j, fabs(tmp), orthogonal_tol) ;
	flag = 0;
      } 
    }

  if ( (verbose_flag == 1) && (flag == 1) )
    log_file << "passed" << endl ;

  return flag ;
}

void qr_decomp(vector<vector<double> > & grad, vector<vector<double> > & t_vec)
{
  int lda , lwork, info ;

  if (verbose_flag == 1)
    log_file << "QR decomposition..." ;

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
      t_vec[i][j] = qr_grad_mat[(k+i)*d + j] ;
  }

  if (verbose_flag == 1)
    log_file << "finished" << endl ;

  check_orthognality(grad, t_vec) ;
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
  int flag ;
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
  flag = 1 ;

  // Newton iteration method 
  while (eps > eps_tol)
  {
    if (verbose_flag == 1)
      log_file << "\t Step " << step << ",\teps = " << eps << endl ;

    if (step + 1 > newton_max_step) 
    {
      flag = - newton_max_step ;
      break ;
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
	mat_a[j*k+i] = s ;
      }

    // right hand side
    for (int i = 0 ; i < k; i ++)
      linear_sol[i] = -tmp_rhs[i] ;

    nrhs = 1; lda = k ; ldb = k ; 
    dgesv_(&k, &nrhs, mat_a, &lda, ipiv, linear_sol, &ldb, &info) ;

    if (info != 0)
    {
      printf("Warning: return value of the linear solver DGESV is wrong: info=%d\n", info) ;
      flag = - (step + 1);
      break ;
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

  if (flag == 1)
  {
    // Newton iteration is successful if we are here
    x0 = tmp_x ;

    if (verbose_flag == 1)
      log_file << "\t Step " << step << ",\teps = " << eps << endl ;

    flag = step ;

    // only count when the convergence is achieved 
    tot_success_newton_step += step ;
    newton_solver_converge_counter ++ ;
  }

  free(vec_a) ;

  tot_newton_step += step ;
  // count how many times we call this function
  newton_solver_counter ++ ;

  return flag ;
}

void compute_v_prime(vector<double> & state, vector<double> & tmp_state, vector<vector<double> > & t_vec) 
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
    coeff[i] = vec_dot(state_diff, t_vec[i]) ;

  for (int i =0; i < d; i ++)
  {
    s = 0 ;
    for (int j = 0 ; j < d-k; j ++)
      s+= coeff[j] * t_vec[j][i];
    v_vec_prime[i] = s;
  }
}

int main ( int argc, char * argv[] ) 
{
  char buf[50] ;
  ofstream out_file ;
  int idx ;
  int newton_success_flag, determinant_flag ;

  double accept_prob , tmp , v_norm_1, v_norm_2 ;

  clock_t start , end ;

  read_config() ;

  mcmc_flag = 1 ;
  orthogonal_tol = 1e-14 ;
  determinant_tol = 1e-5 ;

  size_s = sqrt(2.0 * h_mcmc) ;

  // tolerance of reversibility check
  reverse_tol = 1e-8 ;

  // total newton steps performed so far
  tot_newton_step = 0 ; 

  // total newton steps performed so far, only when converged
  tot_success_newton_step = 0 ;

  newton_solver_counter = 0 ;
  newton_solver_converge_counter = 0 ;

  // rejection counters 
  forward_newton_counter = 0 ;
  backward_newton_counter = 0 ;
  metropolis_counter = 0 ;
  reverse_check_counter = 0 ;
  determinant_counter = 0 ;

  mean_xi_distance = 0 ;

  init_rand_generator();

  printf("SO(N), N=");
  cin >> N ;
  d = N * N ;
  k = N * (N+1) / 2 ;

  n_v_bins = 500; 
  len_v_bound = 5.0 * size_s * sqrt(d-k) ;
  len_v_bin_width = len_v_bound / n_v_bins ;
  bin_len_v_newton_converged.resize(n_v_bins) ;
  bin_len_v_newton_not_converged.resize(n_v_bins) ;

  min_v_norm_not_converged = 1e8 ;
  max_v_norm_converged = -1.0 ;

  //start from the identity matrix
  state.resize(d, 0) ;
  for (int i = 0 ; i < N; i++)
    state[i * N + i] = 1.0 ;

  tmp_state.resize(d) ;
  y_state.resize(d) ;

  // get memory for several global variables
  allocate_mem() ;

  sprintf(buf, "./data/log_ex4_mcmc_son_%d.txt", N) ;
  log_file.open(buf) ;

  // count the time
  start = clock() ;

  printf("\nSO(%d),\td=%d\tk=%d\n", N, d, k) ;
  printf("n=%d\t output_step =%d \n", n, output_every_step ) ;
  printf("h (mcmc) =%.3e,\t size_s=%.3e\n", h_mcmc, size_s) ;

  // for the initial state, compute the Jaccobi (gradient) matrix of xi at current state 
  grad_xi(state, grad_vec) ;

  // for the initial state, fine the orthogonal vectors of tangent space by QR decomposition 
  qr_decomp(grad_vec, tangent_vec_array) ;

  int progress ;
  progress = 1 ;

  for (int i = 0 ; i < n ; i ++)
  {
    if (verbose_flag == 1)
      log_file << "\n==== Generate " << i << "th sample...\n" ;

    if (i >= progress * 0.01 * n)
    {
      end = clock() ;
      printf("\ni=%d, n=%d, %.1f\%% finished.\n", i, n, progress * 1.0 ) ;
      printf("  %.1fsec, %.1esec per step,  remaining time: %4.1fmin\n  Newton success rate=%.2f,  average xi=%.2e\n", (end-start) * 1.0 / CLOCKS_PER_SEC, (end - start) * 1.0 / CLOCKS_PER_SEC / i, 
	  (end - start) * 1.0 / CLOCKS_PER_SEC * (n-i) / (i * 60.0), (i-forward_newton_counter) * 1.0 / i, mean_xi_distance / i );
      progress ++ ;
    }

    trace_series[i] = trace(state) ;

    // randomly generate a vector on the (d-k)-dimensional tangent space
    generate_v(v_vec) ;

    // move along tangent vector v
    for (int ii = 0 ; ii < d; ii ++)
      tmp_state[ii] = state[ii] + v_vec[ii] ;

    mean_xi_distance += distance_to_levelset(tmp_state) ;

    if (verbose_flag == 1)
      log_file 
	<< "1. Move by v : distance = " 
	<< distance_to_levelset(tmp_state) 
	<< endl 
	<< "2. Forward Newton starts:\n" ;

    // projection by newton method
    newton_success_flag = projection_by_Newton(grad_vec, tmp_state) ;

    tmp = sqrt( vec_dot(v_vec, v_vec) )  ;
    idx = int (tmp / len_v_bin_width) ;
    if (idx > n_v_bins-1) idx = n_v_bins - 1;

    if (newton_success_flag < 0) // increase the counter, if we didn't find a new state
    {
      forward_newton_counter ++ ;

      if (verbose_flag == 1)
	log_file << "Rejected by Newton projection, step = " << -newton_success_flag << endl ;

      bin_len_v_newton_not_converged[idx] ++ ;

      if (tmp < min_v_norm_not_converged) min_v_norm_not_converged = tmp ;

      continue;
    }

    if (tmp > max_v_norm_converged) max_v_norm_converged = tmp ;

    bin_len_v_newton_converged[idx] ++ ;

    // check the determinant of the new state
    determinant_flag = check_determinant(tmp_state) ;

    if (determinant_flag == 0)
    {
      determinant_counter ++ ;

      if (verbose_flag == 1)
	log_file << "Rejected due to negative determinant. " << endl ;

      continue ;
    }

    /* 
     * we found a NEW state: tmp_state. Now we do Metropolis-Hasting 
     *
     */

    if (verbose_flag == 1)
      log_file << "\tNewton steps = " << newton_success_flag << endl << "3. Distance after projection = " << distance_to_levelset(tmp_state) << endl ;

    y_state = tmp_state ;

    // compute the Jaccobi (gradient) matrix of xi at the NEW state 
    grad_xi(y_state, grad_vec_y) ;

    // fine the orthogonal vectors of tangent space (NEW state) by QR decomposition 
    qr_decomp(grad_vec_y, tangent_vec_array_y) ;

    // decide vector v'
    compute_v_prime(state, y_state, tangent_vec_array_y) ;

    v_norm_1 = vec_dot(v_vec, v_vec);
    v_norm_2 = vec_dot(v_vec_prime, v_vec_prime);

    // density ratio of two vectors v and v'
    accept_prob = exp(-(v_norm_2 - v_norm_1) * 0.5 / (size_s * size_s) ) ;

    /* 
     * In fact, in this special case, we can show that |v|^2=|v'|^2.
     *  
     */
    if (verbose_flag == 1)
    {
      log_file << "4. |v'|^2=" << v_norm_2 << "\t|v|^2=" << v_norm_1 << "\taccept_prob =" << accept_prob << endl ;

      if ( fabs(v_norm_1 - v_norm_2) > 1e-6)
	log_file << "Warning: the norms of v' and v should be the same!\n" ;
    }

    if (accept_prob > 1.0) accept_prob = 1.0 ;

    tmp = ranf() ;
    if (tmp > accept_prob) // rejected 
    {
      metropolis_counter ++ ;
      continue ;
    }

    /*
     *
     * We have passed the Metropolis-Hasting check. 
     * In the following, we do reversibility check.
     *
     * Again, in this special case, in fact we can show that the reversibility
     * check below is not needed.
     *
     */

    // move state y along tangent vector v'
    for (int ii = 0 ; ii < d; ii ++)
      tmp_state[ii] = y_state[ii] + v_vec_prime[ii] ;

    if (verbose_flag == 1)
      log_file << "5. Backward Newton starts: " << endl ;

    // projection by newton method
    newton_success_flag = projection_by_Newton(grad_vec_y, tmp_state) ;

    if (newton_success_flag < 0) // increase the counter, if we didn't find a new state
    {
      backward_newton_counter ++ ;

      if (verbose_flag == 1)
	log_file << "Rejected by Newton projection, step = " << -newton_success_flag << endl ;

      continue;
    }

    if (verbose_flag == 1)
      log_file << "\tNewton steps = " << newton_success_flag << endl;

    // check whether tmp_state == state .
    tmp = 0;
    for (int ii = 0; ii < d; ii++)
      tmp += fabs(tmp_state[ii] - state[ii]) * fabs(tmp_state[ii] - state[ii]) ;

    if (sqrt(tmp) > reverse_tol)
    {
      reverse_check_counter ++ ;

      if (verbose_flag == 1)
	log_file << "Rejected by reversibility check " << sqrt(tmp) << '>' << reverse_tol << endl ;
    }

    // move to the new state
    state = y_state ;
    // update the gradients and the tangent vectors
    grad_vec = grad_vec_y ;
    tangent_vec_array = tangent_vec_array_y ;
  }

  int tot_rej ;

  printf("\naverage Newton iteration steps = %.2f, ", tot_newton_step * 1.0 / newton_solver_counter );
  if (newton_solver_converge_counter > 0) 
    printf("average Newton steps to reach convergence = %.2f\n", tot_success_newton_step * 1.0 / newton_solver_converge_counter ) ;
  else printf("\n") ;

  printf("\naverage xi distance = %.3e\n", mean_xi_distance * 1.0 / n) ;
  tot_rej = forward_newton_counter + backward_newton_counter + reverse_check_counter + metropolis_counter + determinant_counter ;
  printf("\nRejection rate: Forward\tReverse\tReversibility\tMetrolis\tDeterminant\tTotal\n") ;
  printf("\t\t %.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", forward_newton_counter * 1.0 / n, backward_newton_counter * 1.0 / n, reverse_check_counter *1.0/n, metropolis_counter * 1.0 / n, determinant_counter * 1.0 / n, tot_rej * 1.0 / n) ;

  analysis_data_and_output() ;

  printf("\nmax |v|, when Newton converges : %.4e\n", max_v_norm_converged) ; 
  printf("min |v|, when Newton doesn't converge : %.4e\n", min_v_norm_not_converged) ; 


  sprintf(buf, "./data/ex4_mcmc_v_norm_counter_%d.txt", N) ;
  out_file.open(buf) ;

  out_file << n - forward_newton_counter << ' ' << forward_newton_counter << ' ' << size_s << ' ' << len_v_bound << ' ' << n_v_bins << endl ;

  for (int i = 0 ; i < n_v_bins ; i++)
  {
    out_file << bin_len_v_newton_converged[i] << ' ';
  }
  out_file << endl ;

  for (int i = 0 ; i < n_v_bins ; i++)
  {
    out_file << bin_len_v_newton_not_converged[i] << ' ';
  }
  out_file << endl ;

  out_file.close() ;
  log_file.close() ;

  end = clock() ;
  printf("\n\nRuntime : %4.2f sec.\n\n", (end - start) * 1.0 / CLOCKS_PER_SEC ) ;

  deallocate_mem() ;


  return 0; 
}
