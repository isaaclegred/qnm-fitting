#!/usr/bin/env python
import SetupData
import SetupTrial
from scipy.optimize import minimize, least_squares
import qnm
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
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
                            resolution_level=6,sampling_routine=None, num_samples=None,
                            include_noise=False, plot_confidence_intervals=False,
                            plot_waveforms = True):
    try:
        Yl2m2 = SetupData.get_Yl2m2(data_dir + "/Lev" + str(resolution_level) + \
                                    "/rhOverM_Asymptotic_GeometricUnits_CoM.h5")
    except:
        print("It's possible this file does not exist, try setting up a data directory by using G\
        etAndSetupSXSData.sh")
    start_and_end_frame = SetupData.get_frames_from_offset_and_steps(Yl2m2, offset, num_steps)
    start_frame = start_and_end_frame[0]
    end_frame = start_and_end_frame[1]
    max_frame = start_frame + offset
    ksc = qnm.cached.KerrSeqCache(init_schw=False)
    # The the number of parameters, 2 times the number of modes + Nonlinear params
    numparams = 2*num_modes + 2
    # The initial dimensionless spin of the black hole to use
    # to generate the modes, will be minimized
    A = .75
    M = 1
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
    signal = np.stack((subset_Yl2m2[:, 2],
                       subset_Yl2m2[:, 1]))
    if (include_noise):
        noise  = SetupData.approximate_noise_of_subset(data_dir,
                                               offset, num_steps, 10, included_points)
    else:
        noise = np.ones((len(included_points)))
    # Define residuals (the quantites which when squared and summed give the cost)
    def Residuals(x, noise,  params0, params2):
        target = params0
        grid = params2/x[15]
        trial = SetupTrial.construct_trial_from_grid(x, grid)
        residuals = np.concatenate([target[0, :] - trial[0, :], target[1, :] - trial[1, :]])
        flattened_noise = np.concatenate([noise[:, 0], noise[:, 1]])
        return residuals/flattened_noise
    x0 = np.ones(numparams)
    x0[14] = A
    x0[15] = M
    # Actual fitting
    print("Fitting starting at time " + str(Yl2m2[start_frame, 0]) + " M")
    print("Peak Strain is around 3696.37 M ")
    print("Fitting until " + str(Yl2m2[end_frame,0]) +" M")
    # Compute the best fit given the cost function
    X = least_squares(Residuals, x0 , args=(noise, signal, start_grid), ftol=20**-14, gtol = 10**-15)
    if (plot_waveforms):
        # Plot the Fitted waveform versus the waveform predicted by Numerical
        # Relativity
        plt.plot(start_grid/X['x'][15], signal[1], '-g', label = "NR")
        plt.plot(start_grid/X['x'][15], (SetupTrial.construct_trial_from_grid(X['x'],
                                            start_grid/X['x'][15]))[1], label = "Fit")
        plt.xlabel(r"time $(M)$")
        plt.ylabel(r"$h_{+}$")
        #plt.yscale("log")
        plt.legend()
        plt.savefig("FitMassAndSpin.png")
        plt.figure()
    if(plot_confidence_intervals):
        J = np.matrix(X['jac'])
        H = np.transpose(J)*J
        Avals = np.linspace(-.01, .01, 1000)
        Mvals = np.linspace(-.01, .01, 1000)
        result = np.zeros((len(Avals), len(Mvals)))
        for i in range(len(Avals)):
            for j in range(len(Mvals)):
                result[i,j] = H[14,14]*Avals[i]**2 + 2*H[14,15]*Avals[i]*Mvals[j] + \
                    H[15,15]*Mvals[j]**2 + X['cost']
        print(X['cost'])
        print(len(included_points))
        cs = plt.contour(Mvals, Avals, log(result)/log(10),levels=10)
        plt.ylabel(r"$a-a_{best\,\, fit}$")
        plt.xlabel(r"$M - M_{best\,\, fit}$")
        plt.title("Taylor Expansion of Cost function near minimum")
        plt.clabel(cs)
        plt.clabel(cs)
        #plt.axhline(y=0, color='r', linestyle='-')
        #plt.axvline(x=0, color='r', linestyle='-')
        fig = plt.gcf()
        fig.set_size_inches(10,10)
        plt.savefig("ConfidenceIntervals.png")
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
    return parser.parse_args()
if __name__ == "__main__":
    input_args = global_parse_args()

    fit_qnm_modes_to_signal(input_args.data_dir, input_args.offset, input_args.num_steps,
                            input_args.num_modes, input_args.resolution_level, input_args.sampling_routine,
                            input_args.num_samples, input_args.include_noise,
                            input_args.plot_confidence_intervals,
                            input_args.plot_waveforms)
