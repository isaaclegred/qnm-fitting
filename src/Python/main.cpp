//#include "VectorHelpers.hpp"
#include <iostream>
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

#include "nr3devel.h"
//#include "sep_marquardt.hpp"
namespace bp = boost::python;
namespace np = boost::python::numpy;
// // Get the shape of a numpy array
std::vector<size_t> get_shape(const np::ndarray& A);
//   std::vector<size_t> Shape;
//   for(size_t i = 0; i < A.get_nd(); i++)
//     Shape.push_back(A.shape(i));
//   return Shape;
// }
// // Get the size (as in the length of the A.data()) of an numpy array
size_t get_size(const np::ndarray& A);
//   auto shape = get_shape(A);
//   return  std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<>() );
// }
// // Get a vectorf from a numpy array
VecDoub get_vector_from_np(np::ndarray& A);
//   VecDoub v(get_size(A));
//   auto data = A.get_data();
//   for(size_t j = 0; j < v.size(); ++j){
//     v[j] = data[j];
//   }
//   return v;
// }

int main(){
  int arr_size = 5;
  std::vector<double> v{1,2,3,4,5};
  std::cout << "woeking" << "\n";
  auto data =  v.data();
  std::cout << "working 48" << "\n";
  np::dtype dt = np::dtype::get_builtin<double>();
  bp::tuple shape = bp::make_tuple(arr_size);
  bp::tuple stride = bp::make_tuple(sizeof(double));
  bp::object own;
  std::cout << "working 52" << "\n";
  np::ndarray data_ex1 = np::from_data(data,dt, shape,stride,own);
  std::cout << "working 53" << "\n";
  auto dv = get_vector_from_np(data_ex1);
  std::cout << "working 54" << "\n";
  for (auto element : dv){
    std::cout << element << "\n";
  }
}
