import numpy as np
import qnm

# Using
sin = np.sin
cos = np.cos
real = np.real
imag = np.imag
exp = np.exp

def get_test_funcs(A, grid):
    """
    Construct the analytic QNM's of a black hole with dimless spin A
    """
    ksc = qnm.cached.KerrSeqCache(init_schw=False)
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
def construct_trial(x, test_func):
    """
    Given fitting parameters `x` and known QNM's `test_func` return the
    trial waveform which will be used to compare to the signal data.
    """
    trial = np.zeros((2,len(test_func[0])) )
    for i in range(7):
        trial += np.stack((x[2*i]*test_func[2*i] - x[2*i+1]*test_func[2*i+1],
        x[2*i+1]*test_func[2*i] + x[2*i]*test_func[2*i+1]))
    return trial

def get_flattened_test_funcs(A, grid):
    """
    Construct the analytic QNM's of a black hole with dimless spin A
    """
    ksc = qnm.cached.KerrSeqCache(init_schw=False)
    test_funcs = []
    for i in range(7):
        mode_seq = ksc(s = -2, l = 2, m = 2, n = i)
        freq = mode_seq(a = A)[0]
        hplus = sin(real(freq)*grid)*exp(imag(freq)*grid)
        htimes = cos(real(freq)*grid)*exp(imag(freq)*grid)
        test_funcs.append(np.concatenate([hplus, htimes]))
    return test_funcs

def construct_flattened_trial(x, test_func):
    """
    Construct the trial function in the form of [h_+, h_x] as a contiguous
    array.
    """
    trial = np.zeros((2,len(test_func[0])) )
    for i in range(7):
        trial += np.concatenate([x[2*i]*test_func[2*i] - x[2*i+1]*test_func[2*i+1],
                                 x[2*i+1]*test_func[2*i] + x[2*i]*test_func[2*i+1]])
    return trial
