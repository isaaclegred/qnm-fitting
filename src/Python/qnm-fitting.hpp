#include <vector>
#include "nr3devel.h"
// This is a particular example of a fitting function to be used in an instance of
// the Sep_marquardt class.
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
  std::vector<double> funs(2*num_modes);
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
  std::cout << "funs at t = " << t << "are" <<  funs[0] << " " << funs[1];
  return funs;
}
// This is a bit complicated, because we have real and imaginary parts, we will have 2 times as many
// data points as time values, so we need to treat them as two parts of the same signal, so we will
// treat the imaginary values as the value of the function at time points passed some
// time marker, the `end_time`
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
    std::cout << "There are " << basis_funs.size() << " basis_funs" << "\n";
    // `jac` should have two rows total for real and imaginary parts
    size_t num_modes = f_vals.size()/2;
    std::cout << num_modes << "\n";
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
      for (size_t k = 0; k < 2; k++){
        for (size_t l = 0; l < 2; l++ ){
          std::cout << jac[k][l];
        }
      }
    }
  }
};
