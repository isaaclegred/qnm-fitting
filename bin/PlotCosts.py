import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import qnm
from scipy.optimize import minimize, least_squares
import SetupData
ksc = qnm.cached.KerrSeqCache(init_schw=False)
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp
from MinimizeGivenMa import *

def plot_minimal_costs(Yl2m2_data, offset, num_steps,  a_bounds, M_bounds, a_steps=30, M_steps=30,
                       target_a=None, target_M=None, num_modes=7, precessing=False, spin=1):
    """
    Plot the "cost" associated with the best fit of the coefficients for
    values of a in a_bounds = (a_min, a_max), and M_bounds = (M_min, M_max)
    using a_steps and M_steps grid_points respectively.
    """
    Avals = np.linspace(a_bounds[0], a_bounds[1], a_steps)
    Mvals = np.linspace(M_bounds[0], M_bounds[1], M_steps)
    result = np.zeros((len(Avals), len(Mvals)))
    for i in range(len(Avals)):
        for j in range(len(Mvals)):
            result[i,j] = best_linear_fit_cost(Yl2m2_data, Avals[i], Mvals[j], spin)
    plt.contourf(Mvals,Avals, np.log(result)/np.log(10), levels=60, extend = "both")
    if(target_a):
        plt.axhline(y=target_a, color='r', linestyle='-')
    if(target_M):
        plt.axvline(x=target_M, color='r', linestyle='-')
    plt.xlabel(r"$M_f (M_{i})$")
    plt.ylabel(r"$\chi_f$")
    plt.savefig("GWMinimumCostsManda.png")

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
        '--lower-a ',
        help="The minimum value of a to plot",
        type=float,
        default=.5,
        dest='lower_a'
    )
    parser.add_argument(
        '--upper-a ',
        help="The maximum value of a to plot",
        type=float,
        default=.95,
        dest='upper_a'
    )
    parser.add_argument(
        '--lower-M ',
        help="The minimum value of M to plot",
        type=float,
        default=.7,
        dest='lower_M'
    )
    parser.add_argument(
        '--upper-M ',
        help="The maximum value of M to plot",
        type=float,
        default=1,
        dest='upper_M'
    )
    parser.add_argument(
        '--a-steps ',
        help="The number of steps in the a direction",
        type=int,
        default=30,
        dest='a_steps'
    )
    parser.add_argument(
        '--M-steps ',
        help="The number of steps in the M direction",
        type=int,
        default=30,
        dest='M_steps'
    )
    parser.add_argument(
        '--resolution-level ',
        help="The refinement level to use",
        type=float,
        default=3,
        dest='resolution_level'
    )
    parser.add_argument(
        '--target-a ',
        help="The target dimensionless spin",
        type=float,
        default=None,
        dest='target_a'
    )
    parser.add_argument(
        '--target-M ',
        help="The target mass",
        type=float,
        default=None,
        dest='target_M'
    )
    parser.add_argument(
        "--precessing",
        help="True if the binary system was precessing/ if the spin is not aligmed with z",
        type=bool,
        dest='precessing',
        default=False)
    return parser.parse_args()
if __name__ == "__main__":
    input_args = global_parse_args()
    data_dir = input_args.data_dir
    resolution_level = input_args.resolution_level
    # We will tolerate both having and not having a `/` at the end of data_dir
    slash = ""
    if data_dir[-1] != '/':
        slash += "/"
    if input_args.precessing:
        print("A precessing waveform is being used")
        Yl2m2  = SetupData.get_corrected_2_2(data_dir + slash + "Lev" + \
                                             str(resolution_level) + \
                                    "/rhOverM_Asymptotic_GeometricUnits_CoM.h5")
    else:
        Yl2m2 = SetupData.get_Yl2m2(data_dir + slash + "Lev" + str(resolution_level) + \
                                    "/rhOverM_Asymptotic_GeometricUnits_CoM.h5")
    a_bounds = (input_args.lower_a, input_args.upper_a)
    M_bounds = (input_args.lower_M, input_args.upper_M)
    plot_minimal_costs(Yl2m2, input_args.offset, input_args.num_steps,  a_bounds, M_bounds, input_args.a_steps,
                    input_args.M_steps, input_args.target_a, input_args.target_M,
                       input_args.num_modes, input_args.precessing)
