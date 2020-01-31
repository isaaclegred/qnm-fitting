#include "fgauss_sep.h"
struct fgauss{
  fgauss(std::vector<double> params){};
  void operator()(const Doub x, VecDoub_I &a, VecDoub &u, MatDoub_O &duda){
      return fgauss_sep(x, a, u, duda);}
    void operator()(const Doub t, VecDoub_I & params, Doub & f_eval, VecDoub_O \
                    &jac){}
  };
using fitting_fun = fgauss;
