import numpy as np
from numpy.linalg import svd
from scipy.optimize import minimize
import h5py
import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt
import SetupData
import SetupTrial
import qnm

# Using
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp

h5_file = "SXSDATA/SXS0305Lev0/rhOverM_Asymptotic_GeometricUnits_CoM.h5"
high_res_file  = "SXSDATA/SXS0305Lev6/rhOverM_Asymptotic_GeometricUnits_CoM.h5"
# Define the needed things
Yl2m2 = SetupData.get_Yl2m2(h5_file)
offset = 20
total_steps = 800
max_strain_step = SetupData.find_maxs(Yl2m2)[0]
start_frame  = max_strain_step - offset
end_frame = start_frame + total_steps
grid = Yl2m2[start_frame : end_frame, 0] - Yl2m2[start_frame, 0]* \
np.ones(Yl2m2[start_frame:end_frame,0].size)
signal = np.concatenate((Yl2m2[start_frame:end_frame,2] , Yl2m2[start_frame:end_frame,1]))
print signal.size
# Estimating error by the differnce there are probably better ways to do this, but
# just as a start to have some measure of the error
sigmas = np.ones(signal.size)*np.std(SetupData.get_diff(SetupData.get_Yl2m2(high_res_file), Yl2m2, offset, total_steps))
def get_design_mat(grid, test_funcs, sigmas):
    # The grid is only the number of time steps long, everything else is doubled
    design_mat = np.zeros((2*len(grid), len(test_funcs)))
    print design_mat.shape
    print sigmas.shape,
    for i, func in enumerate(test_funcs):
        print func.shape
        design_mat[:,i] =  func
    return np.matrix(design_mat)
def get_cost(A, grid, signal, sigmas):
    a = A[0]
    test_funcs = SetupTrial.get_flattened_test_funcs(a, grid)
    A = get_design_mat(grid, test_funcs, sigmas)
    U, W, V_dag = svd(A)
    d = np.transpose(np.matrix(signal/sigmas))
    return (np.linalg.norm((1-U*np.transpose(U))*d))**2
A  = .5
Avals  = []
while A < .8:
    Avals.append([A])
    A += .01
Costs = []
for aval in Avals:
    Costs.append(get_cost(aval,grid, signal, sigmas))
plt.plot(Avals, Costs)
plt.savefig("Costs.png")
A0 = .75
args = (grid, signal, sigmas)
X = minimize(get_cost, A0,  args, bounds = [(0,1)])
