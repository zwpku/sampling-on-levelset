#include "ex4.h"

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

void analysis_data_and_output()
{

  ofstream out_file ;
  double s, dd, sigma, tau ;
  char buf[100];
  int idx ;

  // Step 1, output the time series of the trace
  if (mcmc_flag == 0)
    sprintf(buf, "./data/ex4_no_mcmc_traj_%d.txt", N) ;
  else 
    sprintf(buf, "./data/ex4_mcmc_traj_%d.txt", N) ;

  out_file.open(buf) ;
  out_file << n / output_every_step << endl ;

  for ( int i = 0 ; i < n ; i += output_every_step )
    out_file << trace_series[i] << ' ' ;
  out_file << endl ;

  out_file.close() ;

  // Step 2, output the distributions
  
  // divied [-trace_b, trace_b] to n_bins with equal width
  bin_width = 2.0 * trace_b  / n_bins ;
  counter_of_each_bin.resize(n_bins, 0) ;
  for (int i = 0 ; i < n ; i ++)
  {
    idx = int ((trace_series[i] + trace_b) / bin_width) ;
    counter_of_each_bin[idx] ++ ;
  }
  if (mcmc_flag == 0)
    sprintf(buf, "./data/ex4_no_mcmc_counter_%d.txt", N) ;
  else 
    sprintf(buf, "./data/ex4_mcmc_counter_%d.txt", N) ;

  out_file.open(buf) ;
  out_file << n << ' ' << trace_b << ' ' << n_bins << endl ;
  for (int i = 0 ; i < n_bins ; i++)
    out_file << counter_of_each_bin[i] << ' ';
  out_file << endl ;
  out_file.close() ;

  // Step 3, compute mean, variance, autocorelation...
  
  vector<double> cor ;
  cor.resize(max_lag+1) ;
  // compute the mean value of the data
  s = 0 ;
  for (int i = 0 ; i < n ; i ++)
    s += trace_series[i] ;
  s /= n ;
  // substract the mean value 
  for (int i = 0 ; i < n ; i ++)
    trace_series[i] -= s;

  // correlation for different lag time 
  for (int i = 0 ; i <= max_lag ; i ++)
  {
    cor[i] = 0.0 ;
    for (int j = 0 ; j < n - max_lag; j ++)
      cor[i] += trace_series[j] * trace_series[i + j] ;
    cor[i] /= (n - max_lag) ;
  }

  dd = cor[0] ; 
  for (int i = 1 ; i <= max_lag; i ++)
    dd += 2 * cor[i] ;
  sigma = sqrt(dd / n) ;
  tau = dd / cor[0] ;

  printf("n=%d\tmax_lag=%d\tmean=%.4e\tsigma=%.4e\ttau=%.4e\n", n, max_lag, s, sigma, tau) ;

  if (mcmc_flag == 0)
    sprintf(buf, "./data/ex4_no_mcmc_acf_%d.txt", N) ;
  else 
    sprintf(buf, "./data/ex4_mcmc_acf_%d.txt", N) ;

  out_file.open(buf) ;
  out_file << max_lag << ' ' << s << ' ' << sigma << ' ' << tau << endl ;
  for (int i = 0 ; i <= max_lag ; i ++)
    out_file << cor[i] << ' ' ;
  out_file << endl ;

  out_file.close() ;
}

