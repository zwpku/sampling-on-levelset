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

extern int n, verbose_flag, n_bins, output_every_step, maximal_try_number ;
extern double dt, eps_tol, h, trace_b ;

extern "C" {
  extern void dgetrf_(int *, int *, double *, int *, int *, int *) ;
}



