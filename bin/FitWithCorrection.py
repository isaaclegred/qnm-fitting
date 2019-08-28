#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import qnm
from scipy.optimize import minimize, least_squares
import SetupData
import SetupTrial
ksc = qnm.cached.KerrSeqCache(init_schw=False)
try:
    Yl2m2 = SetupData.get_Yl2m2("/Users/isaaclegred/qnm-fitting/SXSDATA0305/Lev0/rhOverM_Asymptotic_GeometricUnits_CoM.h5")
except:
   print("It's possible this file does not exist, try setting up a data directory by using GetAndSetupSXSData.sh")
plot_waveforms = True
plot_confidence_intervals = True
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp
log = np.log

# Fitting the Mass
# These are the parameters that must be set for the Fitting
# NONLINEAR FITTING SPIN
offset = 10
steps  = 750
start_and_end_frame = \
    SetupData.get_frames_from_offset_and_steps(Yl2m2, offset, steps)
start_frame = start_and_end_frame[0]
end_frame = start_and_end_frame[1]
# The the number of parameters, 2 times the number of modes + Nonlinear params
numparams = 30
# The initial dimensionless spin of the black hole to use
# to generate the modes, will be minimized
A = .75
M = 1

# Get the time grid the problem will be analyzed on, note that part of
# fitting the mass of the black hole means changing the physical values
# of the time at the time points, but topologically the grid stays the same
start_grid = Yl2m2[start_frame : end_frame, 0] - Yl2m2[start_frame, 0]*np.ones(Yl2m2[start_frame:end_frame,0].size)
# Get the target signal for the fitting, in this case it comes pre-flattened, which
# doesn't change much, but just means some care must be taken in plotting.
signal = np.concatenate((Yl2m2[start_frame:end_frame,2] , Yl2m2[start_frame:end_frame,1]))
noise = SetupData.approximate_noise("/Users/isaaclegred/qnm-fitting/SXSDATA0305", offset,steps, 10)
flattened_noise = np.zeros((2*len(start_grid)))
for i in range(start_grid.size):
    for j in range(0,2):
        flattened_noise[2*i+j]  =  noise[i,j]
# Define residuals
def Residuals(x, flat_noise,  signal_waveform, unaltered_grid):
    cost = 0;
    target = signal_waveform
    grid = unaltered_grid/x[15]
    test_func = SetupTrial.get_test_funcs(x[14], grid)
    residuals = target - SetupTrial.construct_flattened_trial(x, test_func)
    return residuals/flat_noise
# Construct an initial guess for the configuration of parameters
x0 = np.ones(numparams)
x0[14] = A
x0[15] = M
# Actual fitting
print("Fitting starting at time " + str(Yl2m2[start_frame, 0]) + " M")
print("Peak Strain is around 3696.37 M ")
print("Fitting until " + str(Yl2m2[end_frame,0]) +" M")
# Compute the best fit given the cost function
X = least_squares(Residuals, x0 , args=(flattened_noise, signal, start_grid), ftol=20**-15, gtol = 10**-15)
if (plot_waveforms):
    # Plot the Fitted waveform versus the waveform predicted by Numerical
    # Relativity
    plt.plot(start_grid/X['x'][15], signal[0: int(len(signal) / 2)], '-g', label = "NR")
    test_funcs = SetupTrial.get_test_funcs(X['x'][14], start_grid/X['x'][15])
    cor_test_funcs = SetupTrial.get_spheroidal_correction_funcs(X['x'][14], start_grid/X['x'][15], l_prime = 3)
    plt.plot(start_grid/X['x'][15], (SetupTrial.construct_flattened_trial(X['x'], test_funcs) + SetupTrial.construct_flattened_trial_correction(X['x'], cor_test_funcs))[0:int(len(signal) / 2)],
    label = "Fit")
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
            result[i,j] = H[14,14]*Avals[i]**2 + 2*H[14,15]*Avals[i]*Mvals[j] + H[15,15]*Mvals[j]**2 + X['cost']
    print(X['cost'])
    print(end_frame - start_frame)
    cs = plt.contour(Mvals, Avals, log(result)/log(10),levels=10)
    plt.ylabel(r"$a-a_{best\,\, fit}$")
    plt.xlabel(r"$M - M_{best\,\, fit}$")
    plt.title("Taylor Expansion of Cost function near minimum")
    plt.clabel(cs)
    #plt.axhline(y=0, color='r', linestyle='-')
    #plt.axvline(x=0, color='r', linestyle='-')
    fig = plt.gcf()
    fig.set_size_inches(10,10)
    plt.savefig("ConfidenceIntervals.png")
