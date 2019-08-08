import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import qnm
from scipy.optimize import minimize, least_squares
import SetupData
ksc = qnm.cached.KerrSeqCache(init_schw=False)
Yl2m2 = SetupData.get_Yl2m2("/Users/isaaclegred/qnm-fitting/SXSDATA0305/Lev6/rh\
OverM_Asymptotic_GeometricUnits_CoM.h5")
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp
from MinimizeGivenMa import *
Yl2m2 = SetupData.get_Yl2m2("/Users/isaaclegred/qnm-fitting/SXSDATA0305/Lev6/rhOverM_Asymptotic_GeometricUnits_CoM.h5")
def plot_minimal_costs(Yl2m2_data, a_bounds, M_bounds, a_steps=30, M_steps=30):
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
    plt.axhline(y=0.691, color='r', linestyle='-')
    plt.axvline(x=0.951, color='r', linestyle='-')
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
    return parser.parse_args()
if __name__ == "__main__":
    input_args = global_parse_args()
    a_bounds = (input_args.lower_a, input_args.upper_a)
    M_bounds = (input_args.lower_M, input_args.upper_M)
    plot_minimal_costs(Yl2m2, a_bounds, M_bounds, input_args.a_steps,
                       input_args.M_steps)
