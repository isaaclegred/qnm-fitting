#!/usr/bin/env python
import SetupData
import SetupTrial
from scipy.optimize import minimize, least_squares
import qnm
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from spherical_functions import LM_index as lm
import quaternion
import scri
from FitWithLigoNoise import fit_qnm_modes_to_signal as fit_modes

# I will fit using the function fit_qnm_modes_to_signal, and will have a 2 parameter fit
# over the possible orientations of the z-axis

def fit_precessing_waveform(data_dir, offset, num_steps, num_modes=7,
                            resolution_level=6, sampling_routine=None, num_samples=None,
                            include_noise=False, plot_confidence_intervals=False,
                            plot_waveforms=True, target_spin=None,
                            target_mass=None, save_name="GW", a_guess=None,
                            M_guess=None):
    # We will tolerate both having and not having a `/` at the end of data_dir
    slash = ""
    if data_dir[-1] != '/':
        slash += "/"
    path =  data_dir + slash + \
            "Lev" + str(resolution_level) + "/rhOverM_Asymptotic_GeometricUnits_CoM.h5"
    h = scri.SpEC.read_from_h5(path+"/Extrapolated_N4.dir")
    current_a = .6
    if a_guess:
        current_a = a_guess
    current_M = .95
    if M_guess:
        current_M = M_guess
    # Define Cost Funciton
    def Cost(Q, current_a, current_M):
        q1 = Q[0]
        q2 = Q[1]
        q = quaternion.quaternion(1,q1,q2,0)
        q = q/np.norm(q)
        rotor = np.concatenate([np.array([q.real]),q.imag ])
        h_rot = h.transform(frame_rotation = rotor)
        l, m = 2, 2
        Yl2m2 = np.transpose(np.stack([h_rot.t, np.real(h_rot.data[:, lm(l, m, h.ell_min)]),
                                          np.imag(h_rot.data[:, lm(l, m, h.ell_min)])]))
        X = fit_modes(data_dir,Yl2m2, offset, num_steps, num_modes, sampling_routine,
                      num_samples,
                  include_noise, False,False,
                  target_spin, target_mass, save_name, current_a, current_M)
        # Update the guesses to reflect the most recent fit data.
        current_a = X["x"][-2]
        current_M = X["x"][-1]
        return X["cost"]
    T0 = np.array([0,0])
    T = minimize(lambda Q :  Cost (Q, current_a, current_M) , T0, tol = num_steps*4000)
    return T
def global_parse_args():
    """
    Parse the command line arguments
    """
    import argparse as ap
    parser = ap.ArgumentParser()
    parser.add_argument(
        "--data-dir ",
        help="The directory to find the data in",
        type=str,
        dest='data_dir'
    )
    parser.add_argument(
        "--offset",
        help="The number of steps before the peak strain to start the fitting",
        type=int,
        dest='offset'

    )
    parser.add_argument(
        "--num-steps",
        help="The number of steps between the beginning and ending of the fitting",
        type=int,
        dest='num_steps'
    )
    parser.add_argument(
         "--num-modes",
         help="The number of of QNM overtones to use in the fitting",
         type=int,
         dest='num_modes',
         default=7
    )
    parser.add_argument(
        "--resolution-level",
        help="The resolution level of the data to be used to fit against  (0 to 6)",
        type=int,
        dest='resolution_level',
        default=3

    )
    parser.add_argument(
        "--sampling-routine",
        help="If not None, then instead of using the whole data consisting of " \
        + "num_steps steps, it will sample" +\
        " `num_samples` samples from the points, using the sampling routine given "\
        + "sampling_routine, possible options include `Every-other` and `Start-heavy`"\
        +" (for picking more values near the start for the analysis)",
        type=str,
        dest='sampling_routine',
        default=None
    )
    parser.add_argument(
        "--num-samples",
        help="The number of samples to draw, if a sampling routine is specified",
        type=int,
        dest='num_samples',
        default = None
    )
    parser.add_argument(
        "--include-noise",
        help="Whether to attempt to esitmate errors using a model of injected LIGO noise",
        type=bool,
        dest='include_noise',
        default=True
    )
    parser.add_argument(
        "--plot-confidence-intervals",
        help="Whether to plot confidence intervals and save a .png of the plot",
        type=bool,
        dest='plot_confidence_intervals',
        default=True
    )
    parser.add_argument(
        "--plot-waveforms",
        help="Whether to plot the fitted and target waveform",
        type=bool,
        dest='plot_waveforms',
        default=True
    )
    parser.add_argument(
        "--target-spin",
        help="Target spin, if it is  known, to be plotted with conf. intervals",
        type=float,
        dest='target_spin',
        default=None
    )
    parser.add_argument(
        "--target-mass",
        help="Target mass, if it is  known, to be plotted with conf. intervals",
        type=float,
        dest='target_mass',
        default=None
    )
    parser.add_argument(
        "--save_name",
        help="The name to prepend to the save files, without extension",
        type=str,
        dest='save_name',
        default="GW"
    ),
    parser.add_argument(
        "--a-guess",
        help="If prior information is available about the spin, a guess can be entered of its magnitude",
        type=float,
        dest='a_guess',
        default=.6
    ),
    parser.add_argument(
        "--M-guess",
        help="A guess for the mass of the remnant black hole",
        type=float,
        dest='M_guess',
        default=.95
    )
    return parser.parse_args()
if __name__ == "__main__":
    input_args = global_parse_args()
    T = fit_precessing_waveform(input_args.data_dir, input_args.offset,
                            input_args.num_steps,
                            input_args.num_modes, input_args.resolution_level,
                            input_args.sampling_routine,
                            input_args.num_samples, input_args.include_noise,
                            input_args.plot_confidence_intervals,
                            input_args.plot_waveforms,
                            target_spin = input_args.target_spin,
                            target_mass=input_args.target_mass,
                            save_name=input_args.save_name)
    print(T)
