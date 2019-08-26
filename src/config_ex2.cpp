#include "ex2.h"
#include <iostream>
#include "libconfig.h++"
 
using namespace libconfig;

template <class T>
int read_value(Config & cfg, string  name, T & val)
{
  try
  {
    cfg.lookupValue(name, val);
  }
  catch(const SettingNotFoundException &nfex)
  {
    cerr << "No " << name << " setting in configuration file." << endl;
    return -1;
  }
  return 0;
}

/*
 *
 * Read parameters from the file: ex2.cfg 
 *
 */

int read_config() 
{
  Config cfg;

  try {
    cfg.readFile("ex2.cfg");
  }
  catch (const FileIOException &fioex)
  {
    cerr << "I/O error while reading file." << endl;
    return -1;
  }
  catch(const ParseException &pex)
  {
    cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
	    << " - " << pex.getError() << endl ;
    return -1;
  }

  if ( read_value(cfg, string("n"), n ) < 0 )
    return -1;

  if ( read_value(cfg, string("n_bins"), n_bins) < 0 )
    return -1;

  // maximal step for Newton method
  if ( read_value(cfg, string("newton_max_step"), newton_max_step) < 0 )
    return -1;

  if ( read_value(cfg, string("h"), h) < 0 )
    return -1;

  if ( read_value(cfg, string("h_mcmc"), h_mcmc) < 0 )
    return -1;

  if ( read_value(cfg, string("verbose_flag"), verbose_flag) < 0 )
    return -1;

  if ( read_value(cfg, string("trace_b"), trace_b) < 0 )
    return -1;

  if ( read_value(cfg, string("max_lag"), max_lag) < 0 )
    return -1;

  if ( read_value(cfg, string("output_every_step"), output_every_step) < 0 )
    return -1;

  // convergence check for Newton or ODE solver
  if ( read_value(cfg, string("eps_tol"), eps_tol) < 0 )
    return -1;

  // for SDE-based method
  if ( read_value(cfg, string("maximal_try_number"), maximal_try_number) < 0 )
    return -1;

  // for SDE-based method, step-size in the ODE solver
  if ( read_value(cfg, string("dt"), dt) < 0 )
    return -1;

  return 0 ;
}

