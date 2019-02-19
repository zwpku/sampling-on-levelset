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

  extern void dgetrf_(int *, int *, double *, int *, int *, int *) ;
}

const double pi = atan(1) * 4 ;

extern int n, N, n_bins , newton_max_step, verbose_flag , output_every_step ;

extern double trace_b, h ;


