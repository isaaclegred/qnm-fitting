import numpy as np
import qnm

# Using
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp

def get_test_funcs(A, grid, nummodes=7, spin=1):
    """
    Construct the analytic QNM's of a black hole with dimless spin A, evaluated at the
    times contained in grid.
    """
    ksc = qnm.cached.KerrSeqCache(init_schw=False)
    test_funcs = []
    for i in range(nummodes):
        for j in range(2):
            mode_seq = ksc(s = -2, l = 2, m = spin*2, n = i)
            freq = mode_seq(a = A)[0]
            if j:
                test_funcs.append(sin(real(freq)*grid)*exp(imag(freq)*grid))
            else :
                test_funcs.append(cos(real(freq)*grid)*exp(imag(freq)*grid))
    return test_funcs
def construct_trial(x, test_func, nummodes=7, spin=1):
    """
    Given fitting parameters `x` and known QNM's `test_func` return the
    trial waveform which will be used to compare to the signal data.
    """
    trial = np.zeros((2,len(test_func[0])) )
    for i in range(nummodes):
        trial += np.stack((x[2*i]*test_func[2*i] - x[2*i+1]*test_func[2*i+1],
        x[2*i+1]*test_func[2*i] + x[2*i]*test_func[2*i+1]))
    return trial
def construct_trial_from_grid(x, grid, nummodes = 7, spin = 1):
    """
    Given fitting parameters `x` and known QNM's `test_func` return the
    trial waveform which will be used to compare to the signal data.
    """
    A = x[2*nummodes]
    test_func = get_test_funcs(A, grid, nummodes, spin)
    trial = np.zeros((2,len(test_func[0])) )
    for i in range(nummodes):
        trial += np.stack((x[2*i]*test_func[2*i] - x[2*i+1]*test_func[2*i+1],
                           x[2*i+1]*test_func[2*i] + x[2*i]*test_func[2*i+1]))
    return trial
def construct_complex_trial_from_grid(x, grid, nummodes=7):
    """
    Given fitting parameters `x` and known QNM's `test_func` return the
    trial waveform which will be used to compare to the signal data.
    """
    A = x[2*nummodes]
    test_func = get_test_funcs(A, grid, nummodes)
    trial = np.zeros((len(test_func[0])), dtype=np.complex128 )
    for i in range(nummodes):
        trial += x[2*i]*test_func[2*i] - x[2*i+1]*test_func[2*i+1] + 1j*(
                           x[2*i+1]*test_func[2*i] + x[2*i]*test_func[2*i+1])
    return trial

def construct_flattened_trial(x, test_func, nummodes=7):
    """
    Construct the trial function in the form of [h_+, h_x] as a contiguous
    array.
    """
    trial = np.zeros((2*len(test_func[0])))
    for i in range(nummodes):
        trial += np.concatenate([x[2*i]*test_func[2*i] - x[2*i+1]*test_func[2*i+1],
                                 x[2*i+1]*test_func[2*i] + x[2*i]*test_func[2*i+1]])
    return trial

def get_spheroidal_correction_funcs(A, grid, l_prime = 3, nummodes=7):
    """
    The fact that the QNMs are computed using spheroidal harmonics means that in order
    to compare to spherical harmonic data, we need to include contributions from
    other spheroidal harmonics that mix into the l=2, m=2 spherical harmonic
    """
    ksc = qnm.cached.KerrSeqCache(init_schw=False)
    test_funcs = []
    for i in range(nummodes):
        for j in range(2):
            mode_seq = ksc(s = -2, l = l_prime, m = 2, n = i)
            freq = mode_seq(a = A)[0]
            # Compute the mixing, see qnm documentation for more details
            C = mode_seq(a = A)[2]
            if j:
                test_funcs.append(real(C[0])*sin(real(freq)*grid)*exp(imag(freq)*grid) - \
                imag(C[0])*cos(real(freq)*grid)*exp(imag(freq)*grid))
            else :
                test_funcs.append(real(C[0])*cos(real(freq)*grid)*exp(imag(freq)*grid) + \
                imag(C[0])*sin(real(freq)*grid)*exp(imag(freq)*grid))

    return test_funcs
def construct_trial_correction(x, cor_test_funcs, nummodes=7):
    """
    Construct the trial function in the form of [h_+, h_x] as a contiguous
    array.  `x` needs to have 2*7 + 2 + 2*7 = 30 entries
    """
    trial = np.zeros((2, len(cor_test_funcs)) )
    for i in range(nummodes):
        trial += np.stack([x[2*i + 16]*test_func[2*i+16] - \
                           x[2*i+17]*test_func[2*i+17],
                           x[2*i + 17]*test_func[2*i + 16] + \
                           x[2*i + 16]*test_func[2*i+17]])
    return trial
def construct_flattened_trial_correction(x, cor_test_funcs):
    """
    Construct the trial function in the form of [h_+, h_x] as a contiguous
    array.  `x` needs to have 2*7 + 2 + 2*7 = 30 entries
    """
    trial = np.zeros((2*len(cor_test_funcs[0])) )
    for i in range(7):
        trial += np.concatenate([x[2*i + 16]*cor_test_funcs[2*i] - \
                                 x[2*i+17]*cor_test_funcs[2*i+1],
                                 x[2*i + 17]*cor_test_funcs[2*i] + \
                                 x[2*i + 16]*cor_test_funcs[2*i+1]])
    return trial
