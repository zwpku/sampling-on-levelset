#include "mcmc.h"

// total steps
int n ;
int N, d, k ;
int tot_newton_step, tot_success_newton_step ;
int n_bins ;

double eps_tol, reverse_tol, size_s ;

double mean_xi_distance ;

int newton_max_step ;

int forward_newton_counter, metropolis_counter, backward_newton_counter, reverse_check_counter, determinant_counter ;

double bin_width, trace_b ;

vector<double> state, v_vec, v_vec_prime, y_state, tmp_state, counter_of_each_bin ;

vector<vector<double> > grad_vec , tangent_vec_array, grad_vec_y, tangent_vec_array_y ;

// arrays for QR decomposition
double * qr_tau, * qr_work , * qr_grad_mat ;

double * mat_a, * linear_sol, * mat_x ; 
int *ipiv ;

// used for debug
double orthogonal_tol, determinant_tol ;

int verbose_flag ;
ofstream log_file ;

int read_config() ;
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
      log_file << "determinant check failed!\t det=" << s << "\n" ;
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

double vec_dot(vector<double> & v1, vector<double> & v2) 
{
  double s ;
  s = 0 ;
  for (int i = 0 ; i < v1.size(); i ++)
    s += v1[i] * v2[i] ;
  return s;
}

double distance_to_levelset(vector<double> & x)
{
  vector<double> tmp_vec ;
  double eps ;
  tmp_vec.resize(k) ;
  xi(x, tmp_vec) ;

  eps = sqrt(vec_dot(tmp_vec, tmp_vec)) ;
  return eps ;
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
  }

  free(vec_a) ;

  tot_newton_step += step ;

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

  orthogonal_tol = 1e-14 ;
  determinant_tol = 1e-10 ;
  // residual for Newton method 
  eps_tol = 1e-12 ;

  // tolerance of reversibility check
  reverse_tol = 1e-10 ;

  // total newton steps performed so far
  tot_newton_step = 0 ; 

  // total newton steps performed so far, only when converged
  tot_success_newton_step = 0 ;

  // rejection counters 
  forward_newton_counter = 0 ;
  backward_newton_counter = 0 ;
  metropolis_counter = 0 ;
  reverse_check_counter = 0 ;
  determinant_counter = 0 ;

  mean_xi_distance = 0 ;

  // statistics for output 
  // divied [-trace_b, trace_b] to n_bins with equal width
  bin_width = 2.0 * trace_b  / n_bins ;
  counter_of_each_bin.resize(n_bins, 0) ;

  init_rand_generator();

  printf("SO(N), N=");
  cin >> N ;
  d = N * N ;
  k = N * (N+1) / 2 ;

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
  printf("n=%d\n", n ) ;

  // for the initial state, compute the Jaccobi (gradient) matrix of xi at current state 
  grad_xi(state, grad_vec) ;

  // for the initial state, fine the orthogonal vectors of tangent space by QR decomposition 
  qr_decomp(grad_vec, tangent_vec_array) ;

  for (int i = 0 ; i < n ; i ++)
  {
    if (verbose_flag == 1)
      log_file << "\n==== Generate " << i << "th sample...\n" ;

    tmp = trace(state) ;
    idx = int ((tmp + trace_b) / bin_width) ;
    counter_of_each_bin[idx] ++ ;

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

    if (newton_success_flag < 0) // increase the counter, if we didn't find a new state
    {
      forward_newton_counter ++ ;

      if (verbose_flag == 1)
	log_file << "Rejected by Newton projection, step = " << -newton_success_flag << endl ;

      continue;
    }

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
	log_file << "Rejected by reversibility check\n" ;
    }

    // move to the new state
    state = y_state ;
    // update the gradients and the tangent vectors
    grad_vec = grad_vec_y ;
    tangent_vec_array = tangent_vec_array_y ;
  }

  int tot_rej ;

  printf("\naverage Newton iteration steps = %.2f, average Newton steps to reach convergence = %.2f\n", tot_newton_step * 1.0 / n, tot_success_newton_step * 1.0 / (n-forward_newton_counter - backward_newton_counter) ) ;
  printf("\naverage xi distance = %.3e\n", mean_xi_distance * 1.0 / n) ;
  tot_rej = forward_newton_counter + backward_newton_counter + reverse_check_counter + metropolis_counter + determinant_counter ;
  printf("\nRejection rate: Forward\tReverse\tReversibility\tMetrolis\tDeterminant\tTotal\n") ;
  printf("\t\t %.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", forward_newton_counter * 1.0 / n, backward_newton_counter * 1.0 / n, reverse_check_counter *1.0/n, metropolis_counter * 1.0 / n, determinant_counter * 1.0 / n, tot_rej * 1.0 / n) ;

  sprintf(buf, "./data/ex4_mcmc_counter_%d.txt", N) ;
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
