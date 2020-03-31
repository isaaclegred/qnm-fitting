import numpy as np
import matplotlib.pyplot as plt
import scri
import quaternion
from spherical_functions import LM_index as lm # for easily accessing mode data
import spherical_functions
from mpl_toolkits import mplot3d
from matplotlib import animation
import mpl_toolkits.mplot3d.axes3d as p3

def get_waveform(data_dir):
    """
    Get the scri waveform object from a certain data_directory
    """
    path = data_dir + "/rhOverM_Asymptotic_GeometricUnits_COM.h5"
    h = scri.SpEC.read_from_h5(path+"/Extrapolated_N4.dir")
    return h
def animate_coprecessing_frame(h, save_name="prec_axis",  start_frame=1000,
                               end_frame=1000):
    """
    Make an animation of the orientation of the coprecessing frame as
    given by the output of LLDominantEigenvector.
    """
    sp_eigenvector = h.LLDominantEigenvector()[start_frame:-end_frame]
    def create_lines(times):
        dims = 3
        lineData = np.empty((dims, len(times)))
        line_data = []
        return lineData
    def update_lines(num, dataLines, lines):
        for line, data in zip(lines, dataLines):
            line.set_data(data[0:2, 0:num])
            line.set_3d_properties(data[2, 0 : num])
        return lines
    fig = plt.figure()
    fig.set_size_inches(10, 10)
    ax = p3.Axes3D(fig)
    data = [np.transpose(sp_eigenvector)]
    lines = [ax.plot(dat[0, 0:1], dat[1, 0:1],
                    dat[2, 0:1]) [0] for dat in data]
    ax.set_xlim3d([-1, 1])
    ax.set_xlabel('X')
    ax.set_ylim3d([-1, 1])
    ax.set_ylabel('Y')
    ax.set_zlim3d([-1, 1])
    ax.set_zlabel('Z')
    ax.set_title('Precession Z-axis ')

    # Creating the Animation object
    line_ani = animation.FuncAnimation(fig, update_lines,
                                       12*np.arange(len(sp_eigenvector)//12),
                                       fargs=(data, lines),
                                       interval=1000, blit=False)
    line_ani.save(save_name + '.mp4', fps=50)
    plt.show()
