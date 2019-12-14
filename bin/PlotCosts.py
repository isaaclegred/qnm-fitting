import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import qnm
from scipy.optimize import minimize, least_squares
import SetupData
ksc = qnm.cached.KerrSeqCache(init_schw=False)
Yl2m2 = SetupData.get_Yl2m2("/home/isaaclegred/qnm-fitting/GetData/SXS_BBH_2140/Lev3/rh\
OverM_Asymptotic_GeometricUnits_CoM.h5")
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp
from MinimizeGivenMa import *

def plot_minimal_costs(Yl2m2_data, a_bounds, M_bounds, a_steps=30, M_steps=30,
                       target_a=None, target_M=None):
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
            result[i,j] = best_linear_fit_cost(Yl2m2_data, Avals[i], Mvals[j])
    plt.contourf(Mvals,Avals, np.log(result)/np.log(10), levels=60, extend = "both")
    plt.axhline(y=target_a, color='r', linestyle='-')
    plt.axvline(x=target_M, color='r', linestyle='-')
    plt.xlabel(r"$M_f (M_{i})$")
    plt.ylabel(r"$\chi_f$")
    plt.savefig("MinimumCostsManda.png")

def global_parse_args():
    """
    Parse the command line arguments
    """
    import argparse as ap
    parser = ap.ArgumentParser()
    parser.add_argument(
        '--dir ',
        help="Where to find GW data",
        type=str,
        default="/home/isaaclegred/qnm-fitting/GetData/SXS_BBH_0305",
        dest='dir'
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
        default=1.2,
        dest='upper_M'
    )
    parser.add_argument(
        '--a-steps ',
        help="The number of steps in the a direction",
        type=float,
        default=20,
        dest='a_steps'
    )
    parser.add_argument(
        '--M-steps ',
        help="The number of steps in the M direction",
        type=float,
        default=20,
        dest='M_steps'
    )
    parser.add_argument(
        '--ref-level ',
        help="The refinement level to use",
        type=float,
        default=3,
        dest='ref_level'
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
    return parser.parse_args()
if __name__ == "__main__":
    input_args = global_parse_args()
    Yl2m2 = SetupData.get_Yl2m2(input_args.dir + "/Lev"+ str(input_args.ref_level) + \
                                "/rhOverM_Asymptotic_GeometricUnits_CoM.h5")
    a_bounds = (input_args.lower_a, input_args.upper_a)
    M_bounds = (input_args.lower_M, input_args.upper_M)
    plot_minimal_costs(Yl2m2, a_bounds, M_bounds, input_args.a_steps,
                    input_args.M_steps, input_args.target_a, input_args.target_M)
