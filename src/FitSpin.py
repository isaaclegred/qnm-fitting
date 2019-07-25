#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import qnm
from scipy.optimize import minimize, least_squares
from fitting import SetupData
ksc = qnm.cached.KerrSeqCache(init_schw=False)
Yl2m2 = SetupData.get_Yl2m2("SXSDATA/SXS0305Lev2/rhOverM_Asymptotic_GeometricUn\
its_CoM.h5")
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp
# These are the parameters that must be set for the Fitting
# Fitting the Spin as well as the coefficients
# Very Slow and not accurate
start_frame = 12383
end_frame = 13200
# The the number of parameters, 2 times the number of modes
numparams = 15
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp

# Need to use the constraint on a

# Get the time grid the problem will be analyzed on
grid = Yl2m2[start_frame : end_frame, 0] - Yl2m2[start_frame, 0]*np.ones(Yl2m2[start_frame:end_frame,0].size)
# Get the target signal for the fitting
signal = np.stack((Yl2m2[start_frame:end_frame,2] , Yl2m2[start_frame:end_frame,1]))
# Construct a list of test functions to be used for the fitting
def test_funcs(a):
    if a> 1:
        a=.999
    test_funcs = []
    for i in range(7):
        for j in range(2):
            mode_seq = ksc(s = -2, l = 2, m = 2, n = i)
            freq = mode_seq(a)[0]
            if j:
                test_funcs.append(sin(real(freq)*grid)*exp(imag(freq)*grid))
            else :
                test_funcs.append(cos(real(freq)*grid)*exp(imag(freq)*grid))
    return test_funcs
# Define the cost function for the fitting
# Parameters are ordered target, test function generator, grid
def Cost(x, params0 , params1, params2):
    cost = 0;
    target = params0
    test_functions = params1(x[14])
    trial = construct_trial(x, test_functions)
    for i in range(params2.size):
        cost += np.linalg.norm(target[:,i]- trial[:,i])**2
    return cost
# Define the way the parameters to be fit are used to construct a trial solution to compare
# to the analytic solution.
def construct_trial(x, test_functions):
    test_func = test_functions
    trial = np.zeros((2,len(test_func[0])) )
    for i in range(7):
        trial += np.stack((x[i]*test_func[2*i] - x[2*i+1]*test_func[2*i+1], x[2*i+1]*test_func[2*i] + x[2*i]*test_func[2*i+1]) )
    return trial
# Construct an initial guess for the configuration of parameters
x0 = np.ones(numparams)
# The guess for a can't be 1 for computational reasons
x0[numparams-1] = .68
# Set bounds
Bounds = [(None, None)]*14
Bounds.append((0,.9999))
print("Fitting starting at time " + str(Yl2m2[start_frame, 0]) + " M")
print("Peak Strain is around 3696.37 M ")
print("Fitting until " + str(Yl2m2[end_frame,0]) +" M")
# Compute the best fit given the cost function
X = minimize(Cost, x0 ,bounds = Bounds, args=(signal, test_funcs, grid), gtol = 10**-10)
print(X)
# Plot the Fitted waveform versus the waveform predicted by Numerical Relativity
plt.plot(grid, signal[0], '-g', label = "NR")
plt.plot(grid, (construct_trial(X['x'], test_funcs(X['x'][14])))[0], label = "Fit")
plt.xlabel(r"time $(M)$")
plt.ylabel(r"$h_{+}$")
plt.savefig("FitSpin.png")
#plt.yscale("log")
plt.legend()
