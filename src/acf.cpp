#include "acf.h"

int max_lag, n, N, mcmc_flag ;
vector<double> data ;
vector<double> cor ;

int read_config() ;

void load_data_from_file()
{
  ifstream in_file ;
  char buf[100];

  if (mcmc_flag == 0)
    sprintf(buf, "./data/ex4_no_mcmc_traj_%d.txt", N) ;
  else 
    sprintf(buf, "./data/ex4_mcmc_traj_%d.txt", N) ;

  in_file.open(buf) ;
  in_file >> n ;
  data.resize(n) ;
  for (int i = 0 ; i < n; i ++)
    in_file >> data[i] ;
  in_file.close() ;
}

void acf()
{
  double s, d, sigma, tau ;
  char buf[100];

  cor.resize(max_lag+1) ;
  // compute the mean value of the data
  s = 0 ;
  for (int i = 0 ; i < n ; i ++)
    s += data[i] ;
  s /= n ;
  // substract the mean value 
  for (int i = 0 ; i < n ; i ++)
    data[i] -= s;

  // correlation for different lag time 
  for (int i = 0 ; i <= max_lag ; i ++)
  {
    cor[i] = 0.0 ;
    for (int j = 0 ; j < n - max_lag; j ++)
      cor[i] += data[j] * data[i + j] ;
    cor[i] /= (n - max_lag) ;
  }

  d = cor[0] ; 
  for (int i = 1 ; i <= max_lag; i ++)
    d += 2 * cor[i] ;
  sigma = sqrt(d / n) ;
  tau = d / cor[0] ;

  printf("n=%d\tmax_lag=%d\tmean=%.4e\tsigma=%.4e\ttau=%.4e\n", n, max_lag, s, sigma, tau) ;

  ofstream out_file ;

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

int main() 
{
  read_config() ;

  load_data_from_file() ;

  acf() ;

  return 0;
}

