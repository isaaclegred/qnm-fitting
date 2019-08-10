import h5py
import numpy as np
import glob
def get_Yl2m2(file_name):
    '''
    Get the Y_l2_m2 dataset for a particular h5 file
    '''
    h5file = h5py.File(file_name, 'r')
    YLMs = h5file['Extrapolated_N2.dir']
    Yl2m2 = YLMs['Y_l2_m2.dat']
    Yl2m2_data = np.asarray(Yl2m2)
    h5file.close()
    return Yl2m2_data
def get_levels(directory_name, levels = list(range(7)) ):
    """
    Given a data directory of the form produced by `GetAndSetupSXSDATA.sh`
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
    for i in range(2):
        for j in range(len(max_diff[:,i])):
            if j > avg_over:
                error_estimate[j,i] = np.std(max_diff[j - avg_over : j , i])
            else:
                error_estimate[j,i] = np.std(max_diff[j : j + avg_over, i])
    return error_estimate
def get_frames_from_offset_and_steps(strain_data, offset, steps):
    max_step, h_res_max  =  find_maxs(strain_data)
    return (max_step - offset, max_step - offset + steps)
