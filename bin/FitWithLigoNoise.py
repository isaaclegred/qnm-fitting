#!/usr/bin/env python
import SetupData
import SetupTrial
from scipy.optimize import minimize, least_squares
import qnm
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt


# Some useful definitons to shorten the code
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp
log = np.log
# This is an executable designed to be used to fit QNM to gravitational wave signals,
# see the options below to see what options are available for the fitting.
def fit_qnm_modes_to_signal(data_dir, Yl2m2, offset, num_steps, num_modes=7,
                            sampling_routine=None, num_samples=None,
                            include_noise=False, plot_confidence_intervals=False,
                            plot_waveforms=True, target_spin=None,
                            target_mass=None, save_name="GW", a_guess=None
                            M_guess=None):

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
        noise = 100*SetupData.ligo_noise_stacked(data_dir, offset, num_steps, included_points,
                                   True, 65)
        print ("signal shape is", signal.shape)
        print (" noise shape is",  noise.shape)
        signal  += noise
    else:
        noise = np.ones((len(included_points)))
    # Fitting Procedure
    # Define residuals (the quantites which when squared and summed give the cost)
    def Residuals(x, noise,  params0, params2):
        target = params0
        grid = params2/x[numparams - 1]
        trial = SetupTrial.construct_trial_from_grid(x, grid, num_modes, spin=1)
        residuals = np.concatenate([target[0, :] - trial[0, :], target[1, :] - trial[1, :]])
        flattened_noise = np.concatenate([noise[0, :], noise[1, :]])
        return residuals/flattened_noise
    # Guesses for the individual parameters, we need to provide intelligent guesses for the mass and spin
    x0 = np.ones(numparams)
    x0[numparams - 2] = A
    x0[numparams - 1] = M
    lowerbounds =  [-25]*(numparams -2)
    lowerbounds  = lowerbounds + [0,0]
    upperbounds  =  [25]*(numparams -2)
    upperbounds  = upperbounds + [.999, 1]
    # Actual fitting
    print("Fitting starting at time " + str(Yl2m2[start_frame, 0]) + " M")
    print("Peak Strain is around",Yl2m2[max_frame, 0] ,"M")
    print("Fitting until " + str(Yl2m2[end_frame,0]) +" M")
    # Compute the best fit given the cost function

    X = least_squares(Residuals, x0 , args=(noise, signal, start_grid), bounds = (lowerbounds, upperbounds))
    if (plot_waveforms):
        # Plot the Fitted waveform versus the waveform predicted by Numerical
        # Relativity
        plt.plot(start_grid/X['x'][numparams - 1], signal[1], '-g', label = "NR")
        plt.plot(start_grid/X['x'][numparams - 1], (SetupTrial.construct_trial_from_grid(X['x'],
                                                    start_grid/X['x'][numparams -1], num_modes))[1],
                 label = "Fit")
        print("a is fit as ", X['x'][numparams - 2])
        print("M is fit as", X['x'][numparams - 1])
        plt.xlabel(r"time $(M)$")
        plt.ylabel(r"$h_{+}$")
        #plt.yscale("log")
        plt.legend()
        plt.savefig(save_name + "FitMassAndSpin.png")
        plt.figure()
    if(plot_confidence_intervals):
        if (target_spin):
           plt.axhline(X["x"][numparams - 2] - target_spin, "r")
        if(target_mass):
           plt.axvline(X["x"][numparams -1] - target_mass, "r")
        J = np.matrix(X['jac'])
        H = np.transpose(J)*J
        Avals = np.linspace(-.01, .01, 1000)
        Mvals = np.linspace(-.01, .01, 1000)
        result = np.zeros((len(Avals), len(Mvals)))
        for i in range(len(Avals)):
            for j in range(len(Mvals)):
                result[i,j] = H[numparams-2,numparams-2]*Avals[i]**2 + \
                2*H[numparams-2,numparams-1]*Avals[i]*Mvals[j] + \
                    H[numparams-1,numparams-1]*Mvals[j]**2 + X['cost']
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
        plt.savefig(save_name + "ConfidenceIntervals.png")
    return X
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
