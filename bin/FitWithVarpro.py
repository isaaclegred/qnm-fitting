#!/usr/bin/env python
import SetupData
import SetupTrial
import NonGRFitting
import qnm
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import qnm


# Some useful definitons to shorten the code
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp
log = np.log
# This is an executable designed to be used to fit QNM to gravitational wave signals,
# see the options below to see what options are available for the fitting. If `GR_param_guess`
# is set to an array-like pair of a spin and mass estimate, then the predicted QNM frequencies
# will be used to generate the initial guess for the frequencies, if `freq_guess` (array-like,
# twice the length of num_modes) is provided
# then this will override the initial frequency guess from the GR params, and instead substitute
# the values in `freq_guess` as the initial guess. In either case, the variable `perturbation`
# (same type as `freq_guess`) allows the user to perturb the initial guess for the frequencies
# by some (hopefully small amount)
# hold modes is an optional list of modes to hold, only is important if `fit_freqs` is True
 
def fit_qnm_modes_to_signal(data_dir, Yl2m2, offset, num_steps, num_modes=7,
                            sampling_routine=None, num_samples=None,
                            include_noise=True, plot_confidence_intervals=False,
                            plot_waveforms=True, print_results=False, target_spin=None,
                            target_mass=None, save_name="GW", a_guess=None,
                            M_guess=None, spin=1, fit_freqs=False,
                            GR_param_guess=None, perturbation=None,
                            freq_guess = None, hold_modes=None, hold_c_modes=None):

    start_and_end_frame = SetupData.get_frames_from_offset_and_steps(Yl2m2, offset, num_steps)
    start_frame = start_and_end_frame[0]
    end_frame = start_and_end_frame[1]
    max_frame = start_frame + offset
    ksc = qnm.cached.KerrSeqCache(init_schw=False)
    # The the number of parameters, 2 times the number of modes + Nonlinear params
    numparams = 2*num_modes + 2
    # The initial dimensionless spin of the black hole to use
    # to generate the modes, will be minimized
    if a_guess:
        A = a_guess
    else:
        A = .6
    if M_guess:
        M = M_guess
    else:
        M = .95
    # If a sampling routine for the points is specified, get the a subset of the points
    # sampled using the routine, otherwise, use all the points in [start_frame, end_frame)
    if (sampling_routine):
        included_points = SetupData.get_included_points(start_frame, end_frame, max_frame, sampling_routine)
    else:
        included_points = np.arange(start_frame, end_frame)
    subset_Yl2m2 = SetupData.get_data_subset(Yl2m2, included_points)
    # Get the time grid the problem will be analyzed on, note that part of
    # fitting the mass of the black hole means changing the physical values
    # of the time at the time points, but topologically the grid stays the same
    start_grid = subset_Yl2m2[:, 0] - \
        subset_Yl2m2[0, 0]*np.ones(subset_Yl2m2[:, 0].size)
    signal = np.stack((subset_Yl2m2[:, 2],
                       subset_Yl2m2[:, 1]))

    if (include_noise):
        noise = 10*SetupData.ligo_noise_stacked(data_dir, offset, num_steps, included_points,True, 65)
        signal  += noise
    else:
        noise = np.ones((signal.shape))
    
    long_signal = np.concatenate([signal[0, :], signal[1, :]])

    long_noise = np.concatenate([noise[0, :], noise[1, :]])
    end_time = start_grid[-1]
    double_times  = np.concatenate([start_grid, start_grid + end_time])
    parameter_guesses = np.ones((num_modes * 2))
    # We can make an assumption based on GR to get near to the correct paramaters
    def not_none(val):
        return np.any(val) != None
    if not_none(GR_param_guess):
        for index in range(len(parameter_guesses)//2):
            mode_seq = ksc(s = -2, l = 2, m = 2, n = index)
            freq = mode_seq(a = GR_param_guess[0])[0]
            parameter_guesses[2*index] = -imag(freq)/GR_param_guess[1]
            parameter_guesses[2*index + 1] = real(freq)/GR_param_guess[1]
    else:
        for index in range(len(parameter_guesses)//2):
            parameter_guesses[2*index] = .08 + .15*index
            parameter_guesses[2*index + 1] = .55 - .03*index
            
    if not_none(perturbation):
        parameter_guesses += perturbation
    # Actual fitting
    print("Fitting starting at time " + str(Yl2m2[start_frame, 0]) + " M")
    print("Peak Strain is around",Yl2m2[max_frame, 0] ,"M")
    print("Fitting until " + str(Yl2m2[end_frame,0]) +" M")
    # Compute the best fit given the cost function
    fit_qnms = NonGRFitting.get_marquardt(double_times, long_signal, \
                                          np.std(long_noise)*np.ones((len(long_noise))), \
                                          parameter_guesses, \
                                          2*num_modes, np.array([end_time]))
    print("end time is", end_time)
    # Call the fitting method
    if not(fit_freqs):
        for mode in range(num_modes):
            fit_qnms.hold(2*mode, parameter_guesses[2*mode])
            fit_qnms.hold(2*mode + 1, parameter_guesses[2*mode+1])
    # hold only those modes specified
    if(not_none(hold_modes)):
        for mode in hold_modes:
            fit_qnms.hold(2*mode, parameter_guesses[2*mode])
            fit_qnms.hold(2*mode+1, parameter_guesses[2*mode + 1])

    if(not_none(hold_c_modes)):
        for mode in hold_c_modes:
            fit_qnms.holdc(2*mode, parameter_guesses[2*mode])
            fit_qnms.holdc(2*mode+1, parameter_guesses[2*mode + 1])
    fit_qnms.fit()
    # for mode in range(num_modes):
    #     fit_qnms.hold(2*mode, parameter_guesses[2*mode])
    #     fit_qnms.hold(2*mode + 1, parameter_guesses[2*mode+1])
    #     fit_qnms.holdc(2*mode, 0)
    #     fit_qnms.holdc(2*mode  + 1, 0)
    # for i in range(3): 
    #     for fitting_mode in range(num_modes):
    #         fit_qnms.free(2*fitting_mode)
    #         fit_qnms.free(2*fitting_mode + 1)
    #         fit_qnms.freec(2*fitting_mode)
    #         fit_qnms.freec(2*fitting_mode + 1)
    #         for other_mode in range(num_modes):
    #             if (other_mode != fitting_mode):
    #                 fit_qnms.hold(2*other_mode, fit_qnms.get_nl_params()[2*other_mode])
    #                 fit_qnms.hold(2*other_mode + 1, fit_qnms.get_nl_params()[2*other_mode + 1])
    #                 fit_qnms.fit()
    # Access Fitting Data here
    nl = fit_qnms.get_nl_params()
    lin = fit_qnms.get_lin_params()
    cost = fit_qnms.get_chisq()
    fit_qnms.compute_covar()
    flat_C = fit_qnms.get_covar()
    n = int(np.sqrt(len(flat_C)))
    
    C = np.reshape(flat_C, (n,n))
    import copy
    not_returnable_data = {"name": data_dir +" with " + str(num_modes) + " modes" }
    not_returnable_data["offset"] = offset
    not_returnable_data["numsteps"] = num_steps
    not_returnable_data["duration"] = start_grid[-1] - start_grid[0]
    not_returnable_data["linear_params"] = lin
    not_returnable_data["nonlinear_params"] = nl
    not_returnable_data["covariance"] = C
    not_returnable_data["cost"] = cost
    returnable_data = copy.deepcopy(not_returnable_data)
    if print_results:
        print("Decay Rate and Frequency are ", nl)
        print("The sizes of the relative excitations are", lin)
    if (plot_waveforms):
        plt.figure()
        fit  = np.zeros(len(start_grid))
        for j in range(num_modes):
            fit  +=(lin[2*j]*cos(nl[2*j+1]*start_grid) - \
                    lin[2*j+1]*sin(nl[2*j+1]*start_grid))*np.exp(-nl[2*j]*start_grid)
        plt.plot(start_grid, fit, label="Fit")
        plt.plot(start_grid, signal[0, :], label="NR")
        plt.legend()
        plt.show()
    return returnable_data
   

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


    )
    return parser.parse_args()
if __name__ == "__main__":
    # Assume the if the function is being called from the command line there is no plotting utility
    # available.  
    #mpl.use("agg")
    mpl.use("tkagg")
    input_args = global_parse_args()
    data_dir = input_args.data_dir
    resolution_level = input_args.resolution_level
    # We will tolerate both having and not having a `/` at the end of data_dir
    slash = ""
    if data_dir[-1] != '/':
        slash += "/"

    Yl2m2 = SetupData.get_Yl2m2(data_dir + slash + "Lev" + str(resolution_level) + \
                                    "/rhOverM_Asymptotic_GeometricUnits_CoM.h5")
    fit_qnm_modes_to_signal(data_dir, Yl2m2, input_args.offset,
                            input_args.num_steps,
                            input_args.num_modes,
                            input_args.sampling_routine,
                            input_args.num_samples, input_args.include_noise,
                            input_args.plot_confidence_intervals,
                            input_args.plot_waveforms,
                            target_spin = input_args.target_spin,
                            target_mass=input_args.target_mass,
                            save_name=input_args.save_name
    )
