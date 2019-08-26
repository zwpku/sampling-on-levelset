#include "ex3.h"
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
 * Read parameters from the file: ex3.cfg 
 *
 */

int read_config() 
{
  Config cfg;

  try {
    cfg.readFile("ex3.cfg");
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

  if ( read_value(cfg, string("id_mat_a_flag"), id_mat_a_flag ) < 0 )
    return -1;

  if ( read_value(cfg, string("stiff_eps"), stiff_eps) < 0 )
    return -1;

  if ( read_value(cfg, string("n"), n) < 0 )
    return -1;

  // maximal step for Newton method
  if ( read_value(cfg, string("large_step_size_flag"), large_step_size_flag) < 0 )
    return -1;

  return 0 ;
}

