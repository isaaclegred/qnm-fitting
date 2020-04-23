from scipy.optimize import minimize, least_squares
import qnm
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

ksc = qnm.cached.KerrSeqCache(init_schw=False)
def cost_of_a_and_M (a_and_M, found_freqs, num_modes):
    freqs = []
    A = a_and_M[0]
    M  = a_and_M[1]
    for i in range(num_modes):
        mode_seq = ksc(s = -2, l = 2, m = 2, n = i)
        freq = mode_seq(a = A)[0]
        freqs.append(-np.imag(freq))
        freqs.append(np.real(freq))
    freqs = (np.array(freqs)).flatten()
    freqs = freqs/M
    weights = np.array([.001,.001, 1,1])
    return (freqs - found_freqs)/weights


def infer_a_and_m_from_freqs(decay_times_and_freqs, Covariance = None):
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
    X = least_squares(lambda x : cost_of_a_and_M(x, useful_freqs, num_modes), guess, bounds = Bounds, )
    return X
if __name__ == "__main__":
    print(infer_a_and_m_from_freqs(np.array([0.08438951, 0.55736911, 0.21861721, 0.51372357])))
