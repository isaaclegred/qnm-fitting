#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import qnm
from scipy.optimize import minimize, least_squares
import SetupData
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
numparams = 16
# The initial dimensionless spin of the black hole to use
# to generate the modes, will be minimized
A = .75
M = 1

# Get the time grid the problem will be analyzed on, note that part of
# fitting the mass of the black hole means changing the physical values
# of the time at the time points, but topologically the grid stays the same
start_grid = Yl2m2[start_frame : end_frame, 0] - Yl2m2[start_frame, 0]*np.ones(Yl2m2[start_frame:end_frame,0].size)
# Get the target signal for the fitting
signal = np.stack((Yl2m2[start_frame:end_frame,2] , Yl2m2[start_frame:end_frame,1]))
noise = SetupData.approximate_noise("/Users/isaaclegred/qnm-fitting/SXSDATA0305", offset,steps, 10)
# Construct a list of test functions to be used for the fitting
def get_test_funcs(A, grid):
    test_funcs = []
    for i in range(7):
        for j in range(2):
            mode_seq = ksc(s = -2, l = 2, m = 2, n = i)
            freq = mode_seq(a = A)[0]
            if j:
                test_funcs.append(sin(real(freq)*grid)*exp(imag(freq)*grid))
            else :
                test_funcs.append(cos(real(freq)*grid)*exp(imag(freq)*grid))
    return test_funcs
# Define residuals
def Residuals(x, noise,  params0, params2):
    cost = 0;
    target = params0
    grid = params2/x[15]
    trial = construct_trial(x, grid)
    residuals = np.zeros((2*len(grid)))
    flattened_noise = np.zeros((2*len(grid)))
    for i in range(grid.size):
        for j in range(0,2):
            residuals[2*i+j] = (target[j,i]- trial[j,i])
            flattened_noise[2*i+j]  =  noise[i,j]
    return residuals/flattened_noise
# Define the way the parameters to be fit are used to construct a trial
# solution to compare to the analytic solution.
def construct_trial(x, grid):
    A = x[14]
    test_func = get_test_funcs(A, grid)
    trial = np.zeros((2,len(test_func[0])) )
    for i in range(7):
        trial += np.stack((x[2*i]*test_func[2*i] - x[2*i+1]*test_func[2*i+1],
                           x[2*i+1]*test_func[2*i] + x[2*i]*test_func[2*i+1]) )
    return trial
# Construct an initial guess for the configuration of parameters
x0 = np.ones(numparams)
x0[14] = A
x0[15] = M
# Actual fitting
print("Fitting starting at time " + str(Yl2m2[start_frame, 0]) + " M")
print("Peak Strain is around 3696.37 M ")
print("Fitting until " + str(Yl2m2[end_frame,0]) +" M")
# Compute the best fit given the cost function
X = least_squares(Residuals, x0 , args=(noise, signal, start_grid), ftol=20**-15, gtol = 10**-15)
if (plot_waveforms):
    # Plot the Fitted waveform versus the waveform predicted by Numerical
    # Relativity
    plt.plot(start_grid/X['x'][15], signal[1], '-g', label = "NR")
    plt.plot(start_grid/X['x'][15], (construct_trial(X['x'], start_grid/X['x'][15]))[1],
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
