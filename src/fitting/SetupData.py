import h5py
import numpy as np
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
    print(Yl2m2)
    abs_strain = abs(Yl2m2[:,1] + 1j*Yl2m2[:,2])
    max_step = np.where(abs_strain == np.amax(np.max(abs_strain)))
    return max_step[0][0], np.amax((abs_strain))
def get_diff(higher_res, lower_res, offset, steps):
    h_res_max_step, h_res_max  =  find_maxs(higher_res)
    l_res_max_step, l_res_max  =  find_maxs(lower_res)
    h_data = higher_res[1:2, h_res_max_step - offset: h_res_max_step - offset + steps]
    l_data  = lower_res[1:2, l_res_max_step - offset: l_res_max_step - offset + steps]
    return (h_data - l_data)
