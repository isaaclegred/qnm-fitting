#include "nr3devel.h"
#include "numeric"

namespace vec{
  template<typename T>
  void SIZE_CHECK(T a, T b){
    if (a.size() != b.size()){
      throw("Data sizes do not match!");
    }
  }
  double inner_product(VecDoub t, VecDoub s){
    SIZE_CHECK<VecDoub>(t,s);
    double value = 0;
    for(size_t i = 0; i < t.size(); ++i){
      value = t[i] * s[i];
    }
  return value;
  }

  VecDoub& times(double c, VecDoub& t){
    for(size_t i = 0; i < t.size(); ++i){
    t[i] *= c;
    }
    return t;
  }
} // namespace vec
