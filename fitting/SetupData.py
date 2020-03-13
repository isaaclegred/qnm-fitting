import h5py
import numpy as np
import glob
from spherical_functions import LM_index as lm
import quaternion
import scri
# Convert a frequency given in units G=c=M=1 to Hertz, bh_mass is in solar
# masses
def convert_freq_to_hz(bh_mass, freq):
    G = 6.67*10**-11
    c = 3*10**8
    return freq*c**3/G/(bh_mass*1.89*10**30)/(2*np.pi)
# Convert the other way
def convert_from_hz(bh_mass, freq):
    G = 6.67*10**-11
    c = 3*10**8
    return freq/c**3/G/(bh_mass*1.89*10**30)/(2*np.pi)

def get_Yl2m2(file_name, spin_dir = 1):
    '''
    Get the Y_l2_m2 dataset for a particular h5 file
    '''

    h5file = h5py.File(file_name, 'r')
    YLMs = h5file['Extrapolated_N2.dir']
    if spin_dir == 1:
        Yl2m2 = YLMs['Y_l2_m2.dat']
    else:
        Yl2m2 = YLMs['Y_l2_m-2.dat']

    Yl2m2_data = np.asarray(Yl2m2)
    h5file.close()
    return Yl2m2_data
def get_levels(directory_name, levels = list(range(7)) ):
    """
    Given a data directory of the form produced by `SXS`
    and an array like levels containing integer refinement levels, return
    a list of data arrays as would be returned from calling get Yl2m2 on each
    file.
    """
    data_arrays = []
    ref_levels = glob.glob(directory_name + "/Lev*")
    # The levels may not be in order, so we sort them
    ref_levels.sort(key  = lambda name : float(name[-1]))
    for ref_level in ref_levels:
        data_arrays.append(get_Yl2m2(ref_level + "/rhOverM_Asymptotic_GeometricUnits_CoM.h5"))
    return data_arrays
def get_data(file_name, data_name ):
    '''
    Get any particular dataset `data_name` from an h5 file
    '''
    h5file = h5py.File(file_name, 'r')
    YLMs = h5file['Extrapolated_N2.dir']
    Yl2m2 = YLMs[data_name]
    Yl2m2_data = np.asarray(Yl2m2)
    h5file.close()
    return Yl2m2_data

def list_dat_files(file_name):
    '''
    Get a list of all possible datasets in an h5 file
    '''
    h5file = h5py.File(file_name, 'r')
    YLMs = h5file['Extrapolated_N2.dir']
    possible_data = list(YLMs.keys())
    h5file.close()
    return possible_data
def find_maxs(time_data):
    '''
    return a list of times steps where the complex strain magnitude
    is maximized, and a list of strains at these points.
    '''

    Yl2m2 = time_data
    abs_strain = abs(Yl2m2[:,1] + 1j*Yl2m2[:,2])
    max_step = np.where(abs_strain == np.amax(np.max(abs_strain)))
    return max_step[0][0], np.amax((abs_strain))
def get_diff(higher_res, lower_res, offset, steps):
    """
    Return the difference between two signals, one with higher resolution and
    one with lower resolution in order to obtain an estimate on the numerical
    noise in the simulation data.
    """
    # The maximum strain doesn't necessarily appear at the same time step in
    # all resolutions, so we compare the data relative to the max strain
    h_res_max_step, h_res_max  =  find_maxs(higher_res)
    l_res_max_step, l_res_max  =  find_maxs(lower_res)
    h_data = higher_res[ h_res_max_step - offset: h_res_max_step - offset + steps, 1:3]
    l_data  = lower_res[ l_res_max_step - offset: l_res_max_step - offset + steps, 1:3]
    return (h_data - l_data)
def approximate_noise(data_dir, offset, steps, avg_over):
    """
    Return an approximation to the numerical noise of Yl2m2, assigning
    to a timestep the standard deviation of the nearest avg_over steps of
    the result of get_diff
    """
    Yl2m2_arrays = get_levels(data_dir)
    diff_arrays = []
    for resolution_1 in Yl2m2_arrays:
        for resolution_2 in Yl2m2_arrays:
                diff_arrays.append(get_diff(resolution_1, resolution_2, offset, steps))
    max_diff = diff_arrays[0]
    for diff_array in diff_arrays:
        if np.sum(diff_array**2 > np.sum(max_diff**2)):
            max_diff = diff_array
    error_estimate = np.zeros(max_diff.shape)
    print(max_diff.shape)
    for i in range(2):
        for j in range(len(max_diff[:,i])):
            if j > avg_over:
                error_estimate[j,i] = np.std(max_diff[j - avg_over : j , i])
            else:
                error_estimate[j,i] = np.std(max_diff[j : j + avg_over, i])
    return error_estimate
def get_frames_from_offset_and_steps(strain_data, offset, steps):
    """
    Get the frames which are the subset of strain_data with indices
    given by the index of the max strain minus the offset, up to that
    index plus the number of steps, steps
    """
    max_step, h_res_max  =  find_maxs(strain_data)
    return (max_step - offset, max_step - offset + steps)
def get_data_subset(strain_data, included_points):
    """
    Given strain_data, and an array containing the indices of which points
    to include in the final data set, return the desired points
    """

    final_data = np.zeros((len(included_points), 3))
    for j in range(final_data.shape[0]):
        final_data[j, :] = strain_data[included_points[j], :]
    return final_data
def approximate_noise_of_subset(data_dir, offset, steps, avg_over,
                                included_points):
    """
    Return an approximation to the numerical noise of a subset of Yl2m2,
    assigning to a timestep the standard deviation of the nearest avg_over
    steps of the result of get_diff
    """
    full_noise = approximate_noise(data_dir, offset, steps, avg_over)
    noise_subset = np.zeros((len(included_points), 2))
    for index, point in enumerate(included_points):
        noise_subset[index, :]  = full_noise[point - min(included_points),:]
    return noise_subset
def get_every_nth_point(n, start_frame, end_frame):
    """
    Get an array containing every nth integer between integers start_frame and end_frame,
    meant to be used to define variable `included_points` in `FitSubsetWithNoise.py`
    Always include the zeroth point
    """
    return  np.asarray([i for i in range(start_frame,  end_frame) if \
                        (i - start_frame) % n == 0])
def sample_with_distribution(P, numsamples,  start_frame, end_frame):
    """
    Sample without replacement from a distribution of integers in the range
    [start_frame, end_frame] with probabilities (or relative probabilities)
    given by an array P,
    repeat numsamples times and return an array of
    samples, sorted in ascending order
    """
    # Nothing will work if the length of the probability distribution is not
    # the same size as the array of ints being sampled from.
    assert(len(P) == end_frame -start_frame)
    # Normalize probabilites
    P = P/sum(P)
    all_timesteps = np.arange(start_frame, end_frame)
    samples = \
        np.random.choice(all_timesteps, size=(numsamples),
                                     replace=False, p=P )
    samples.sort()
    return samples
def get_included_points(start_frame, end_frame, max_frame, sampling_routine):
    steps = end_frame - start_frame
    if (sampling_routine == "Start-heavy"):
        P = 1/(1+ .01*np.arange(steps))
        return sample_with_distribution(P, int(steps/2), start_frame, end_frame)
    if (sampling_routine == "Every-other"):
        return  (np.arange(start_frame, end_frame))[::2]
#  The distance given, R, is the estimated distance 150914 in units of
#  masses of 150914, it will scale linearly with edistnace and inversely
#  with black hole mass.  This is necessary because everything internally is
#  done in terms of masses of the initial black hole binary "massD"
def ligo_noise_of_subset(data_dir, offset, num_steps, included_points,
                         convert, bh_mass,R=1.24*10**20):
    """
    Use LIGO signal to noise ratios to generate a noise signal designed to
    obscure a GW signal to the point of LIGO design sensitivity

    This model simplifies many points, including the assumption that only noise
    near the ringdown frequencies will be important.  The SNR depends as 1/R
    on the distance from the source to the detector, but numerical data is scaled
    to eliminate this factor, so it is reintroduced here to mimic the effect of
    more distant detections.
    """
    Ligo_Noise = np.loadtxt("ligonoise.dat")
    # Find the range of frequencies which correspond to usual QNM values, and will serve
    # to interfere with the fitting
    critical_freqs = (min(np.where(Ligo_Noise[:, 0] > 240)[0]),
                      max(np.where(Ligo_Noise[:,0]<265)[0]))
    # Take the maximum noise value from the critical range of frequencies
    error_source = max(Ligo_Noise[critical_freqs[0] : critical_freqs[1],
                                  1])*R
    if(convert):
        error_source = convert_from_hz(bh_mass, error_source)
        # Next line is the chain rule term to take into account the fact that
        # the strain error has units sqrt(1/hz), so to convert it to units of
        # 1/M we need to multiply by \sqrt(1/M)/\sqrt(hz), which is the squae
        # root of the thing we just computed
        error_source  = error_source/np.sqrt(error_source)
    noise  = np.random.normal(0, scale=error_source, size=(len(included_points)))
    return noise
def ligo_noise_stacked(data_dir, offset, num_steps, included_points,
                         convert, bh_mass,R=1.24*10**20):
    Ligo_Noise = np.loadtxt(data_dir + "/ligonoise.dat")
    # Find the range of frequencies which correspond to usual QNM values, and will serve
    # to interfere with the fitting
    critical_freqs = (min(np.where(Ligo_Noise[:, 0] > 240)[0]),
                      max(np.where(Ligo_Noise[:,0]<265)[0]))
    # Take the maximum noise value from the critical range of frequencies
    error_source = max(Ligo_Noise[critical_freqs[0] : critical_freqs[1],
                                  1])
    error_source =  error_source**(1/2)*R
    print(error_source)
    noise  = 1/np.sqrt(2)*np.random.normal(0, scale=error_source, size=(2, len(included_points)))
    return noise
def get_l_2_data(file_name):
    all_data = []
    for m in ["-2", "-1", "0", "1", "2"]:
        h5file = h5py.File(file_name, 'r')
        YLMs = h5file['Extrapolated_N2.dir']
        all_data.append(YLMs["Y_l2_m" + m + ".dat"])
    return all_data
def get_corrected_2_2(file_name):
    h = scri.SpEC.read_from_h5(file_name + "/Extrapolated_N2.dir")
    # Get the waveform in the coprecessing frame
    i = quaternion.quaternion(0,1,0,0)
    j = quaternion.quaternion(0,0,1,0)
    k = quaternion.quaternion(0,0,0,1)
    def cross(a,b):
        return 1/2*(a*b - b*a)
    k_hat = quaternion.quaternion(0, 0.25269634778, 0.04202266793, 0.80968828452 )
    k_hat = k_hat/ np.sqrt(k_hat*np.conjugate(k_hat))
    ax = cross(k, k_hat)
    ax = ax/np.sqrt(ax*np.conjugate(ax))
    cos_al = np.real(k*np.conjugate(k_hat))
    alpha = np.arccos(cos_al)
    rotor = np.concatenate([np.array([np.cos(0.5*alpha)]),np.sin(0.5*alpha)*ax.imag])
    h.to_corotating_frame()
    #h_rot = h.transform(frame_rotation=rotor)
    h_rot = h
    return np.transpose(np.stack([h_rot.t,
                                  np.real(h_rot.data[:, lm(2, 2, h_rot.ell_min)]),
                                  np.imag(h_rot.data[:, lm(2, 2, h_rot.ell_min)])]))
