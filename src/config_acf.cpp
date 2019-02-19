#include "acf.h"
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
 * Read parameters from the file: acf.cfg 
 *
 */

int read_config() 
{
  Config cfg;

  try {
    cfg.readFile("acf.cfg");
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

  if ( read_value(cfg, string("N"), N) < 0 )
    return -1;

  if ( read_value(cfg, string("mcmc_flag"), mcmc_flag) < 0 )
    return -1;

  if ( read_value(cfg, string("max_lag"), max_lag) < 0 )
    return -1;

  return 0 ;
}

