#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/stl_iterator.hpp>
#include <cmath>
#include <string>
#include <numeric>
#include <vector>
#include <cstdio>
#include <tuple>
#include <utility>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <iostream>

#include "nr3devel.h"
#include "sep_marquardt.hpp"

namespace bp = boost::python;
namespace np = boost::python::numpy;
// Compute the vector of [Re{\psi_1(t)}, ... Re{\psi_n(t)}, Im{\psi_1(t)}, ...  Im{\psi_n(t)}]
// params are [x_1, ...x_n, y_1, ... y_n, \omega_1, ... \omega_n, \Gamma_1, \Gamma_n], where
// x_i  = A_i cos(\phi_i), where A_1 is the real amplitude, and \phi_1 is the phase
// y_i = A_i sin(\phi_i), this stems from the complex vector with N entries, and converting
// it into a vector of 2N real numbers.  The conversion is c_n, complex amplitude is
// c_n  = A_n e^{i \phi_n} and \Omega_n a complex frequency is \Omega_n =  \omega_n + i \Gamma_n
// c_n e^{i \Omega_n t} = e^{-\Gamma_n t}[x_n cos \omega_n t - y_n sin \omega_n t ,
// x_n sin \omega_n t + y_n cos \omega_n t]
std::vector<double> get_basis_funs(VecDoub params, double t, bool real){
// There are 4 real parameters per basis function, if there are not a whole number of
// basis functions worth of parameters, it should be made known, but
  if (params.size() % 2 !=0){
    printf("Warning! an incorrect  number of parameters have been passed get_basis_funs");
  }
  size_t num_modes = params.size() / 2;
  std::vector<double> funs(num_modes);
  // For each function we get a real and imaginary part
  for (size_t i  = 0; i < 2*num_modes; ++i){
      double Gamma = params[2 * i] ;
      double omega = params[2*i + 1];
      if (real){
        funs[2*i] = exp(-Gamma*t) * cos(omega*t);
        funs[2*i + 1] =  - exp(-Gamma*t) * sin(omega*t);
      } else{
        funs[2*i] = exp(-Gamma*t) * sin(omega*t) ;
        funs[2*i +  1] = exp(-Gamma*t) * cos(omega*t);
      }

  }
  return funs;
}
// This is a bit complicated, because we have real and imaginary parts, we will have 2 times as many
// data points as time values, so we need to treat them as two parts of the same signal, so we will
// treat the imaginary values as the value of the function at
struct fitting_fun{
  double end_time;
  fitting_fun(double time) : end_time(time){};
  // For plugging into the superclass constructor,
  void operator()(const Doub t, VecDoub_I & params, Doub & f_eval, VecDoub_O &jac){}
  void operator()(const Doub  t, VecDoub_I & params , VecDoub_O & f_vals,
           MatDoub_O &  jac){
    // If t  <= end_time then we are computing a real value, so the real parts of the basis
    // functions are used.
    bool real =  not (t > end_time);
    // Put in expressions for the f_vals and jac
    // $$f(t_i)  = \sum_n c_n e^{-i \omega_n t_i}$$
    auto basis_funs =  get_basis_funs(params, t, real);
    // `jac` should have two rows total for real and imaginary parts
    size_t num_modes = f_vals.size()/2;
    for(size_t j = 0; j < num_modes; ++j){
      // The jacobian will be a 2N \times 2N matrix (as there are 2N ``functions",
      // and 2N parameters, the real and imaginary parts of the Omega_j, \omega_j and \Gamma_j),
      // it will only be nonzero on the diagonals of the N \times N sub-block matrices ,
      // The i ``rotates" the real  part into the imaginary part,
      // and the imaginary part into the negative real part

      // If real then these imaginary times -1, if not real, then these are the real functions
      auto other_funs_for_t = get_basis_funs(params, t, not real);
      // derivative of e^{-Gamma_j t} \cos \omega_j t if real or e^{-Gamma_j t} \sin \omega_j t
      // with respect to \Gamma_j
      jac[2*j][2*j] = -t * basis_funs[2*j];
      // derivative of e^{-Gamma_j t} \cos \omega_j t if real or e^{-Gamma_j t} \sin \omega_j t
      // with respect to \omega_j
      jac[2*j][2*j + 1] = -t * basis_funs[2*j + j] ;
      // // derivative of e^{-Gamma_j t} -\sin \omega_j t if real or e^{-Gamma_j t} \cos \omega_j t
      // with respect to \omega_j
      jac[2*j + 1][2*j] = t * basis_funs[2*j + 1 ];
      jac[2*j + 1][2*j + 1] = -t * basis_funs[2*j] ;
    }
  }
};
using varpro = Sep_marquardt<fitting_fun>;

// Get the shape of a numpy array
std::vector<size_t> get_shape(const np::ndarray& A){
  std::vector<size_t> Shape;
  for(size_t i = 0; i < A.get_nd(); i++)
    Shape.push_back(A.shape(i));
  return Shape;
}
// Get the size (as in the length of the A.data()) of an numpy array
size_t get_size(const np::ndarray& A){
  auto shape = get_shape(A);
  return  std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<>() );
}
// Get a vectorf from a numpy array
VecDoub get_vector_from_np(np::ndarray& A){
  VecDoub v(get_size(A));
  auto data = A.get_data();
  for(size_t j = 0; j < v.size(); ++j){
    v[j] = data[j];
  }
  return v;
}


varpro get_marquardt(np::ndarray& xxa, np::ndarray& yya,
                    np::ndarray& ssiga, np::ndarray& aaa,
                    fitting_fun f,
                    const int pp,  const double TOL=1.e-3) noexcept{
    return varpro(get_vector_from_np(xxa),
                         get_vector_from_np(yya),
                         get_vector_from_np(ssiga),
                         get_vector_from_np(aaa),
                         f, pp);
}


BOOST_PYTHON_MODULE(mylib){
  bp::def("get_marquardt", +[](double time, np::ndarray& xxa, np::ndarray& yya,
                               np::ndarray& ssiga, np::ndarray& aaa, const int pp){
                              return get_marquardt(xxa, yya, ssiga, aaa, fitting_fun(time), pp);
                            });
  bp::def("get_numpy", +[](np::ndarray& input){
                          int arr_size = get_size(input);
                          auto v = get_vector_from_np(input);
                          auto data = v.data();
                          np::dtype dt = np::dtype::get_builtin<double>();
                          bp::tuple shape = bp::make_tuple(arr_size);
                          bp::tuple stride = bp::make_tuple(sizeof(double));
                          bp::object own;
                          np::ndarray data_ex1 = np::from_data(data,dt, shape,stride,own);
                          return data_ex1;
                        }
    );
  bp::class_<varpro>("Varpro", bp::no_init)
    .def("holdc", +[](varpro& M, const int i,
                      const double val) {M.holdc(i, val);})
    .def("freec",+[](varpro& M, const int i) {M.freec(i);})
    .def("hold", +[](varpro& M, const int i,
                     const double val) {M.hold(i, val);})
    .def("free", +[](varpro& M, const int i) {M.free(i);})
    .def("compute_covar", +[](varpro& M) { M.compute_covar();})
    .def("fit", +[](varpro& M) { M.fit();});
  bp::def("greet", +[](){
                      std::string s =  "hello";
                      const char* c = s.c_str();
                      printf("%s", c);
                      return bp::str(c) ;
                    });
  Py_Initialize();
  np::initialize();
}
