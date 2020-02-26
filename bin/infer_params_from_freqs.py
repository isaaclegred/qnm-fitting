from scipy.optimize import minimize, least_squares
import qnm
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

ksc = qnm.cached.KerrSeqCache(init_schw=False)
def fit (a_and_M, found_freqs, num_modes):
    freqs = []
    A = a_and_M[0]
    M  = a_and_M[1]
    for i in range(num_modes):
        mode_seq = ksc(s = -2, l = 2, m = 2, n = i)
        freq = mode_seq(a = A)[0]
        freqs.append(3*(1-1/(i+1))**2*np.real(freq))
        freqs.append(-3 *(1-1/(i+1))**2*np.imag(freq))
    (np.array(freqs)).flatten()
    return freqs - found_freqs


def infer_a_and_m_from_freqs(decay_times_and_freqs):
    num_modes = int(len(decay_times_and_freqs)/2)
    sortable_freqs = []
    for i in range(num_modes):
        sortable_freqs.append([decay_times_and_freqs[i], decay_times_and_freqs[i + 1]])
    sortable_freqs.sort(key = lambda x : -x[1])
    useful_freqs = np.ndarray.flatten(np.array(sortable_freqs))
    a_guess = .7
    M_guess = .95
    guess = [a_guess, M_guess]
    Bounds = ([0,0], [.99999,1])
    X = least_squares(lambda x : fit(x, useful_freqs, num_modes), guess, bounds = Bounds, )
    return X
if __name__ == "__main__":
    print(infer_a_and_m_from_freqs(np.array([0.23524088, 0.64771748, 0.18619304, 0.42950507, 0.07782983, 0.55375245])))
