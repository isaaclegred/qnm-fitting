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
// Change the below include to choose a fitting function
#include "qnm-fitting.hpp"
//#include "fgauss.hpp"
namespace bp = boost::python;
namespace np = boost::python::numpy;

using varpro = Sep_marquardt<fitting_fun>;
// using varpro =  Sep_marquardt<fgauss>;
/// Get the shape of a numpy array
std::vector<size_t> get_shape(const np::ndarray& A){
  std::vector<size_t> Shape;
  for(size_t i = 0; i < A.get_nd(); i++)
    Shape.push_back(A.shape(i));
  return Shape;
}
/// Get the size (as in the length of the A.data()) of an numpy array
size_t get_size(const np::ndarray& A){
  auto shape = get_shape(A);
  return  std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<>() );
}
/// Get a std::vector<double> from a numpy array
VecDoub get_vector_from_np(np::ndarray& A){
  VecDoub v(get_size(A));
  double* data = reinterpret_cast<double*>(A.get_data());
  for(size_t j = 0; j < v.size(); ++j){
    v[j] = data[j];
  }
  return v;
}
/// Get a numpy array from a std::vector<double>
np::ndarray get_np_from_vector(const std::vector<double>& v){
  const double* data = v.data();
  np::dtype dt = np::dtype::get_builtin<double>();
  bp::tuple shape = bp::make_tuple(v.size());
  bp::tuple stride = bp::make_tuple(sizeof(double));
  np::ndarray my_zeros = np::zeros(shape, dt);
  double* zero_data =
    reinterpret_cast<double*>(my_zeros.get_data());
  for (size_t i = 0; i < v.size(); ++i) {
    zero_data[i] = data[i];
    std::cout << "data is trying to be np-ified" << data[i] << "\n";
  }
  return my_zeros;
}
varpro get_marquardt(np::ndarray& xxa, np::ndarray& yya,
                    np::ndarray& ssiga, np::ndarray& aaa,
                    fitting_fun f,
                    const int pp,  const double TOL=1.e-1) noexcept{
  std::cout << "there are " << get_size(aaa) <<" nonlinear params " <<"\n";
    return varpro(get_vector_from_np(xxa),
                         get_vector_from_np(yya),
                         get_vector_from_np(ssiga),
                         get_vector_from_np(aaa),
                         f, pp);

 }
// varpro get_marquardt_gauss(np::ndarray& xxa, np::ndarray& yya,
//                     np::ndarray& ssiga, np::ndarray& aaa,
//                     fgauss f,
//                     const int pp,  const double TOL=1.e-1) noexcept{
//   std::cout << "there are " << get_size(aaa) <<" nonlinear params " <<"\n";
//     return varpro(get_vector_from_np(xxa),
//                          get_vector_from_np(yya),
//                          get_vector_from_np(ssiga),
//                          get_vector_from_np(aaa),
//                   f, pp);}

BOOST_PYTHON_MODULE(mylib){
  /// get_marquardt should be implemented in an external .hpp file, and fitting_fun should be an
  /// an alias for the name of struct which contains a callable function as specified in
  /// `Sep_marquardt.hpp`
  bp::def("get_marquardt", +[](np::ndarray& xxa, np::ndarray& yya,
                               np::ndarray& ssiga, np::ndarray& aaa, const int pp,
                                   np::ndarray& params, const double TOL =1e-3){
                              auto times = get_vector_from_np(xxa);
                              auto v = get_vector_from_np(xxa);
                              std::cout << "gettting a marquardt" << "\n";
                              return get_marquardt(xxa, yya, ssiga, aaa,
                                                   fitting_fun(get_vector_from_np(params)),
                                                   pp, TOL);
                               });
  // bp::def("get_marquardt_gauss", +[](np::ndarray& xxa, np::ndarray& yya,
  //                              np::ndarray& ssiga, np::ndarray& aaa, const int pp){
  //                             auto times = get_vector_from_np(xxa);
  //                             auto v = get_vector_from_np(xxa);
  //                             for(auto element : v){std::cout << element << "\n";}
  //                             return get_marquardt_gauss(xxa, yya, ssiga, aaa,
  //                                                  fgauss(), pp);
  //                            });
  bp::def("get_numpy", +[](np::ndarray& input) -> np::ndarray {
                          int arr_size = get_size(input);
                          std::cout << arr_size << '\n';
                          const auto v = get_vector_from_np(input);
                          const double* data = v.data();
                          np::dtype dt = np::dtype::get_builtin<double>();
                          bp::tuple shape = bp::make_tuple(arr_size);
                          bp::tuple stride = bp::make_tuple(sizeof(double));
                          np::ndarray my_zeros = np::zeros(shape, dt);
                          double* zero_data =
                            reinterpret_cast<double*>(my_zeros.get_data());
                          for (size_t i = 0; i < arr_size; ++i) {
                            zero_data[i] = data[i];
                          }
                          return my_zeros;
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
    .def("fit", +[](varpro& M) { M.fit();})
    .def("get_times", +[](varpro& M) -> np::ndarray {
                         return get_np_from_vector(M.x);})
    .def("get_values", +[](varpro& M) -> np::ndarray {
                          return get_np_from_vector(M.y);})
    .def("get_nl_params", +[](varpro& M) -> np::ndarray{
                             return get_np_from_vector(M.a);})
    .def("get_lin_params", +[](varpro& M) -> np::ndarray{
                              return get_np_from_vector(M.c);})
    .def("get_chisq", +[](varpro& M){return M.chisq;});
  bp::def("greet", +[](){
                      std::string s =  "hello";
                      const char* c = s.c_str();
                      printf("%s", c);
                      return bp::str(c) ;
                    });
  Py_Initialize();
  np::initialize();
}
