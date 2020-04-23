#!/usr/bin/env python
import SetupData
import SetupTrial
from scipy.optimize import minimize, least_squares
import qnm
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import NonGRFitting
# Some useful definitons to shorten the code
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp
log = np.log

# This is an executable designed to be used to fit QNM to gravitational wave signals,
# see the options below to see what options are available for the fitting.
def fit_qnm_modes_to_signal(data_dir, offset, num_steps, num_modes=7,
                            resolution_level=6,sampling_routine=None, num_samples =None,
                            include_noise=False, plot_confidence_intervals=False,
                            plot_waveforms=True, target_spin=None,
                            target_mass=None):
    try:
        Yl2m2 = SetupData.get_Yl2m2(data_dir + "/Lev" + str(resolution_level) +"/rhOverM_Asymptotic_GeometricUnits_CoM.h5")
    except:
        print(data_dir + "/Lev" + str(resolution_level) +\
"/rhOverM_Asymptotic_GeometricUnits_CoM.h5")
        print("It's possible this file does not exist, try setting up a data directory by using GetAndSetupSXSData.sh")
    start_and_end_frame = SetupData.get_frames_from_offset_and_steps(Yl2m2, offset, num_steps)
    start_frame = start_and_end_frame[0]
    end_frame = start_and_end_frame[1]
    max_frame = start_frame + offset
    ksc = qnm.cached.KerrSeqCache(init_schw=False)
    # The the number of parameters, 2 times the number of modes + Nonlinear params
    numparams = 4*num_modes
    # Get the time grid the problem will be analyzed on, note that part of
    # fitting the mass of the black hole means changing the physical values
    # of the time at the time points, but topologically the grid stays the same
    start_grid = Yl2m2[start_frame : end_frame, 0] - Yl2m2[start_frame, 0]*np.ones(Yl2m2[
        start_frame:end_frame,0].size)
    # If a sampling routine for the points is specified, get the a subset of the points
    # sampled using the routine, otherwise, use all the points in [start_frame, end_frame)
    if (sampling_routine):
        included_points = SetupData.get_included_points(start_frame, end_frame, max_frame, sampling_routine)
    else:
        included_points = np.arange(start_frame, end_frame)
    subset_Yl2m2 = SetupData.get_data_subset(Yl2m2, included_points)
    start_grid = subset_Yl2m2[:, 0] - \
        subset_Yl2m2[0, 0]*np.ones(subset_Yl2m2[:, 0].size)
    signal = np.stack((subset_Yl2m2[:, 2], subset_Yl2m2[:, 1]))

    if (include_noise):
        noise = 1*SetupData.ligo_noise_stacked(data_dir, offset, num_steps, included_points,
                                   True, 65)
        print ("signal shape is", signal.shape)
        print ("noise shape is",  noise.shape)
        signal  += noise
    else:
        noise = np.ones((len(included_points)))
    # The fitting here is done with an underlying C++ implementation of the varpro algorithm
    # Any variables defined on a time grid need to be doubled
    long_signal = np.concatenate([signal[0, :], signal[1, :]])

    plt.show()
    long_noise = np.concatenate([noise[0, :], noise[1, :]])
    end_time = start_grid[-1]
    double_times  = np.concatenate([start_grid, start_grid + end_time])
    parameter_guesses = (.5) * np.ones((num_modes * 2))
    for index, guess in enumerate(parameter_guesses):
        print("index is", index)
        parameter_guesses[index] -= .02 * index
    fit_qnms = mylib.get_marquardt(double_times, long_signal, np.std(long_noise)*np.ones((len(
        long_noise))), parameter_guesses, 2*num_modes, np.array([end_time]))
    print("end time is", end_time)
    # Call the fitting method
    for mode in range(num_modes):
        fit_qnms.hold(2*mode, .5)
        fit_qnms.hold(2*mode + 1, .5)
        fit_qnms.holdc(2*mode, 0)
        fit_qnms.holdc(2*mode  + 1, 0)
    for fitting_mode in range(num_modes):
        fit_qnms.free(2*fitting_mode)
        fit_qnms.free(2*fitting_mode + 1)
        fit_qnms.freec(2*fitting_mode)
        fit_qnms.freec(2*fitting_mode + 1)
        for other_mode in range(num_modes):
            if (other_mode != fitting_mode):
                fit_qnms.hold(2*other_mode, fit_qnms.get_nl_params()[2*other_mode])
                fit_qnms.hold(2*other_mode + 1, fit_qnms.get_nl_params()[2*other_mode + 1])
        fit_qnms.fit()
    # Access Fitting Data here
    nl = fit_qnms.get_nl_params()
    lin = fit_qnms.get_lin_params()
    print("Decay Rate and Frequency are ", nl)
    plt.figure()
    fit  = np.zeros(len(start_grid))
    for j in range(num_modes):
        fit  +=(lin[2*j+1]*cos(nl[2*j+1]*start_grid) - lin[2*j]*sin(nl[2*j+1]*start_grid))*np.exp(-nl[2*j]*start_grid)
    plt.plot(start_grid, fit, label="Fit")
    plt.plot(start_grid, signal[0, :], label="NR")
    plt.legend()
    plt.show()
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
        default=6

    )
    parser.add_argument(
        "--sampling-routine",
        help="If not None, then instead of using the whole data consisting of num_steps steps, it will " +\
        "sample num_samples samples from the points, using the sampling routine given by sampling_routine",
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
        help="Whether to attempt to esitmate errors using numerical noise from the simulation",
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
    return parser.parse_args()
if __name__ == "__main__":
    input_args = global_parse_args()

    fit_qnm_modes_to_signal(input_args.data_dir, input_args.offset, input_args.num_steps,
                            input_args.num_modes, input_args.resolution_level, input_args.sampling_routine,
                            input_args.num_samples, input_args.include_noise,
                            input_args.plot_confidence_intervals,
                            input_args.plot_waveforms, target_spin = input_args.target_spin,
                            target_mass=input_args.target_mass)
