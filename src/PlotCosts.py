def plot_minimal_costs(a_bounds, M_bounds, a_steps=30, M_steps=30):
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
            result[i,j] = best_linear_fit_cost(Avals[i], Mvals[j])
    plt.contourf(Mvals,Avals, log(result)/log(10), levels=60, extend = "both")
    plt.axhline(y=0.691, color='r', linestyle='-')
    plt.axvline(x=0.951, color='r', linestyle='-')
    plt.xlabel(r"$M_f (M_{i})$")
    plt.ylabel(r"$\chi_f$")
    plt.figure(figsize = (30,26))
    plt.savefig("MinimumCostsManda.png")

def parse_args():
    """
    Parse the command line arguments
    """
    import argparse as ap
    parser = ap.ArgumentParser(
    parser.add_argument(
        '--lower-a '
        required=False
        help="The minimum value of a to plot"
        type=float
        default=.5
    )
    parser.add_argument(
        '--upper-a '
        required=False
        help="The maximum value of a to plot"
        type=float
        default=1
    )
    parser.add_argument(
        '--lower-M '
        required=False
        help="The minimum value of M to plot"
        type=float
        default=.7
    )
    parser.add_argument(
        '--upper-M '
        required=False
        help="The maximum value of M to plot"
        type=float
        default=1.2
    )
    parser.add_argument(
        '--a-steps '
        required=False
        help="The number of steps in the a direction"
        type=float
        default=30
    )
    parser.add_argument(
        '--M-steps '
        required=False
        help="The number of steps in the M direction"
        type=float
        default=30
    )
        
if __name__ == "__main__":
    input_args = parse_args()
        a_bounds = (input_args["lower_a"])
    plot_minimal_costs(a_bounds, M_bounds, a_steps=30, M_steps=30):