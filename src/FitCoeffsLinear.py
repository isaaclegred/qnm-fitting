#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import qnm
from scipy.optimize import least_squares
from fitting import SetupData
from fitting import SetupTrial
ksc = qnm.cached.KerrSeqCache(init_schw=False)
try:
    Yl2m2 = SetupData.get_Yl2m2("/home/isaaclegred/qnm-fitting/SXSDATA0305/Lev0/rhOverM_Asymptotic_GeometricUnits_CoM.h5")
except:
    "It's possible this file does not exist, try setting up a data
    by using GetAndSetupSXSData.sh "
# Use scipy.optimize.least_squares to achieve a much faster fit in the linear case
# These are the parameters that must be set for the Fitting
start_frame = 12300
end_frame = 13200
# The the number of parameters, 2 times the number of modes
numparams = 14
# The dimensionless spin of the black hole to use to generate the modes
A = .75

# Using
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp
# Get the time grid the problem will be analyzed on
grid = Yl2m2[start_frame : end_frame, 0] - Yl2m2[start_frame, 0]*np.ones(Yl2m2[start_frame:end_frame,0].size)
# Get the target signal for the fitting
signal = np.stack((Yl2m2[start_frame:end_frame,2] , Yl2m2[start_frame:end_frame,1]))
# Construct a list of test functions to be used for the fitting
test_funcs = SetupTrial.get_test_funcs(A)
# Define residuals for this particular fit
def Residuals(x, params0, params1, params2):
    cost = 0;
    target = params0
    test_functions = params1
    trial = SetupTrial.construct_trial(x, test_functions)
    residuals = np.zeros((2*len(params2)))
    for i in range(params2.size):
        for j in range(0,2):
            residuals[i+j] = target[j,i]- trial[j,i]
    return residuals

# Construct an initial guess for the configuration of parameters
x0 = np.ones(numparams)
print("Fitting starting at time " + str(Yl2m2[start_frame, 0]) + " M")
print("Peak Strain is around 3696.37 M ")
print("Fitting until " + str(Yl2m2[end_frame,0]) +" M")
# Compute the best fit given the cost function
X = least_squares(Residuals, x0 , args=(signal, test_funcs, grid))
# Plot the Fitted waveform versus the waveform predicted by Numerical Relativity
plt.plot(grid, signal[1], '-g', label = "NR")
plt.plot(grid, (SetupTrial.construct_trial(X['x'], test_funcs))[1], label = "Fit")
plt.xlabel(r"time $(M)$")
plt.ylabel(r"$h_{+}$")
#plt.yscale("log")
plt.legend()
print("Saving plot to " + "NumericalvsFit.png")
plt.savefig("NumericalvsFit.png")
print("The objective function has a value " + str(X["cost"]))
