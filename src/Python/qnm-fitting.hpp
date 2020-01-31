#include "nr3devel.h"
#include "vectorhelpers.hpp"
#include <algorithm>
#include <iterator>
#include <vector>
/// This is a particular example of a fitting function to be used in an instance
/// of the Sep_marquardt class. Compute the vector of \f$[Re{\psi_1(t)}, ...
/// Re{\psi_n(t)}, Im{\psi_1(t)}, ...  Im{\psi_n(t)}]\f$ params are \f$[x_1,
/// ...x_n, y_1, ... y_n, \omega_1, ... \omega_n, \Gamma_1, \Gamma_n]\f$, where
/// \f$x_i  = A_i cos(\phi_i)\f$, where \f$A_1\f$ is the real amplitude, and
/// \f$\phi_1\f$ is the phase \f$y_i = A_i sin(\phi_i)\f$, this stems from the
/// complex vector with \f$N\f$ entries, and converting it into a vector of
/// \f$2N\f$ real numbers.  The conversion is \f$c_n\f$, complex amplitude is
/// \f$ c_n  = A_n e^{i \phi_n}\f$ and \f$\Omega_n\f$ a complex frequency is
/// \f$\Omega_n =  \omega_n + i \Gamma_n\f$ Breaking the expression down into a
/// vector of real and imaginary parts \f$c_n e^{i \Omega_n t} \sim e^{-\Gamma_n
/// t}[x_n cos \omega_n t - y_n sin \omega_n t , x_n sin \omega_n t + y_n cos
/// \omega_n t]\f$
std::vector<double> get_basis_funs(VecDoub params, double t, bool real,
                                   double t_e);

std::vector<double> get_basis_funs(VecDoub params, double t, bool real,
                                   double t_e) {
  // There are 4 real parameters per basis function, if there are not a whole
  // number of basis functions worth of parameters, it should be made known, but
  if (params.size() % 2 != 0) {
    printf("Warning! an incorrect  number of parameters have been passed "
           "get_basis_funs");
  }
  size_t num_modes = params.size() / 2;
  std::vector<double> funs(2 * num_modes);
  // For each function we get a real and imaginary part
  for (size_t i = 0; i < num_modes; i++) {
    double Gamma = params[2 * i];
    double omega = params[2 * i + 1];
    if (real) {
      funs[2 * i] = exp(-Gamma * t) * cos(omega * t);
      funs[2 * i + 1] = -exp(-Gamma * t) * sin(omega * t);
    } else {
      funs[2 * i] = exp(-Gamma * (t - t_e)) * sin(omega * (t - t_e));
      funs[2 * i + 1] = exp(-Gamma * (t - t_e)) * cos(omega * (t - t_e));
    }
  }
  return funs;
}
/// This is a bit complicated, because we have real and imaginary parts, we will
/// have 2 times as many data points as time values, so we need to treat them as
/// two parts of the same signal, so we will treat the imaginary values as the
/// value of the function at time points passed some time marker, the `end_time`
// The jacobian will be a \f$2N \times 2N\f$ matrix (as there are \f$2N\f$
/// ``functions", and \f$2N\f$ parameters, the real and imaginary parts of the
/// \f$\Omega_j,\, \omega_j\f$ and \f$\Gamma_j)\f$, it will only be nonzero on
/// the diagonals of the N \times N sub-block matrices , The i ``rotates" the
/// real  part into the imaginary part, and the imaginary part into the
/// negative real part
/// derivative of with respect to \f$\Gamma_j\f$ \f$e^{-\Gamma_j t} \cos
/// \omega_j t = -t e^{-\Gamma_j t} \cos \omega_j t\f$ if real or
/// \f$e^{-\Gamma_j t} sin \omega_j t  = -t e^{-Gamma_j t} sin \omega_j t\f$
/// derivative of with respect to \f$\omega_j\f$,  \f$e^{-\Gamma_j t} \cos
/// \omega_j t = -t e^{-\Gamma_j t} \sin \omega_j t \f$ if real or \f$
/// e^{-\Gamma_j t} \sin \omega_j t = t e^{-\Gamma_j t} \cos \omega_j t \f$ if
/// imaginary
struct fitting_qnm {
  double end_time;
  // The only parameter for this function is the end_time of the real data, data
  // which comes after the `end_time` are data belonging to the imaginary part
  // of the qnm signal
  fitting_qnm(std::vector<double> contains_final_time) {
    end_time = contains_final_time[0];
  }
  // For plugging into the superclass constructor,
  void operator()(const Doub t, VecDoub_I &params, Doub &f_eval,
                  VecDoub_O &jac) {}
  void operator()(const Doub t, VecDoub_I &params, VecDoub_O &f_vals,
                  MatDoub_O &jac) {
    // If t  <= end_time then we are computing a real value, so the real parts
    // of the basis functions are used.
    bool real = not(t > end_time);
    // Put in expressions for the f_vals and jac
    // f(t_i)  = \sum_n c_n e^{-i \omega_n t_i}
    auto basis_funs = get_basis_funs(params, t, real, end_time);
    f_vals.clear();
    std::copy(basis_funs.begin(), basis_funs.end(), std::back_inserter(f_vals));
    // `jac` should have two rows total for real and imaginary parts
    // std::cout << "the function values look like" << "\n";
    // vec::print(f_vals);
    // The time in the imaginary case is not the true time, it needs to be
    // adjusted
    double t_c = t - end_time * (not(real));
    size_t num_modes = f_vals.size() / 2;
    for (size_t j = 0; j < 2 * num_modes; j++)
      for (size_t i = 0; i < 2 * num_modes; i++)
        jac[j][i] = 0.0;
    for (size_t j = 0; j < num_modes; ++j) {
      // Derivatives with respect to Gamma_j,
      jac[2 * j][2 * j] = -t_c * basis_funs[2 * j];
      jac[2 * j + 1][2 * j] = -t_c * basis_funs[2 * j + 1];
      // In this case the derivatives with respect to omega_j really matter, and
      // the real functions pick up an extra sign flip
      jac[2 * j][2 * j + 1] = t_c * basis_funs[2 * j + 1];
      jac[2 * j + 1][2 * j + 1] = -t_c * basis_funs[2 * j];
    }
  }
};
using fitting_fun = fitting_qnm;
