#include "fgauss_sep.h"
struct fgauss{
    void operator()(const Doub x, VecDoub_I &a, VecDoub &u, MatDoub_O &duda){
      return fgauss_sep(x, a, u, duda);}
    void operator()(const Doub t, VecDoub_I & params, Doub & f_eval, VecDoub_O \
                    &jac){}
  };
