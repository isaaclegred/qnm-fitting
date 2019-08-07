import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import qnm
from scipy.optimize import minimize, least_squares
import SetupData
ksc = qnm.cached.KerrSeqCache(init_schw=False)
Yl2m2 = SetupData.get_Yl2m2("SXSDATA/SXS0305Lev2/rhOverM_Asymptotic_GeometricUn\
its_CoM.h5")
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp
def minimize_given_mass_and_spin(strain_data = Yl2m2, A =.75, M = 1):
    """
    Return the best fit  associated with fitting to the `strain_data`
    signal using a = `A` and M = `M`
    """
    start_frame = 12383
    end_frame = 13200
    # The the number of parameters, 2 times the number of modes
    numparams = 14
    # Get the time grid the problem will be analyzed on
    raw_grid = strain_data[start_frame : end_frame, 0] - strain_data[start_frame, 0]*np.ones(strain_data[start_frame:end_frame,0].size)
    this_grid = raw_grid/M
    # Get the target signal for the fitting
    signal = np.stack((strain_data[start_frame:end_frame,1] , strain_data[start_frame:end_frame,2]))
    # Construct a list of test functions to be used for the fitting
    test_funcs = []
    for i in range(7):
        for j in (True, False):
            mode_seq = ksc(s = -2, l = 2, m = 2, n = i)
            freq = mode_seq(a = A)[0]
            if j:
                test_funcs.append(sin(real(freq)*this_grid)*exp(imag(freq)*this_grid))
            else :
                test_funcs.append(cos(real(freq)*this_grid)*exp(imag(freq)*this_grid))
    # Define residuals
    def Residuals(x, params0, params1, params2):
        cost = 0;
        target = params0
        test_functions = params1
        trial = construct_trial(x, test_functions)
        residuals = np.zeros((2*len(params2)))
        for i in range(params2.size):
            for j in range(0,2):
                residuals[2*i+j] = target[j,i]- trial[j,i]
        return residuals
    # Define the way the parameters to be fit are used to construct a trial solution to compare
    # to the analytic solution.
    def construct_trial(x, test_func):
        trial = np.zeros((2,len(test_func[0])) )
        for i in range(7):
            trial += np.stack((x[2*i]*test_func[2*i] - x[2*i+1]*test_func[2*i+1], x[2*i+1]*test_func[2*i] + x[2*i]*test_func[2*i+1]) )
        return trial
    # Construct an initial guess for the configuration of parameters
    x0 = np.ones(numparams)
    X = least_squares(Residuals, x0 , args=(signal, test_funcs, this_grid), gtol = 10**-15)
    return (X, test_funcs)
def best_linear_fit_cost(a, m):
    """
    Return the cost associated with the best fit to spin and mass.
    """
    (X,test_funcs) = minimize_given_mass_and_spin(A = a, M = m)
    return X["cost"]
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
    plt.contourf(Mvals,Avals, log(result)/log(10), levels=60,extend = "both")
    plt.axhline(y=0.691, color='r', linestyle='-')
    plt.axvline(x=0.951, color='r', linestyle='-')
    plt.xlabel(r"$M_f (M_{i})$")
    plt.ylabel(r"$\chi_f$")
    plt.figure(figsize = (30,26))
    plt.savefig("MinimumCostsManda.png")
