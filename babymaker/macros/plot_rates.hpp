// plot_rates: Plots HLT rates as a function of HT, MET
#ifndef H_PLOT_RATES
#define H_PLOT_RATES

#include <vector>
#include "TString.h"

class var_class {
public:
  var_class(TString ivarname, float iminx, float imaxx, TString ititle);
  TString title, varname;
  float minx, maxx;
};

class sample_class {
public:
  sample_class(TString ifiles, TString ilabel, int icolor);
  TString files;
  TString label;
  int color;
};



#endif
