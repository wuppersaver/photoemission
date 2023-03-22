import numpy as np
import matplotlib.pyplot as plt
import time

def make_ome_hist_plot(files:dict,coords,lower, upper):
    data = {}
    coord_formulas = {
        'x' : r'$|OME(x)|$',
        'y' : r'$|OME(y)|$',
        'z' : r'$|OME(z)|$',
        'xy' : r'$\sqrt{|OME(x)|^2 + |OME(y)|^2}$',
        'xz' : r'$\sqrt{|OME(x)|^2 + |OME(z)|^2}$',
        'yz' : r'$\sqrt{|OME(y)|^2 + |OME(z)|^2}$',
        'xyz' : r'$\sqrt{|OME(x)|^2 + |OME(y)|^2 + |OME(z)|^2}$',
    }
    for key in files.keys():
        data[key] = np.genfromtxt(files[key],skip_header=6).flatten()
    fig,ax = plt.subplots(1,1,figsize=[10,5],dpi=200)
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlabel(coord_formulas[coords])
    ax.set_ylabel('count')
    for key in files.keys():
        ax.hist(data[key],bins=np.logspace(np.log10(lower),np.log10(upper), 150),alpha=0.5,label=key)
    #ax.hist(ome_flat[system],bins=np.logspace(np.log10(lower_value),np.log10(upper_value), 100),alpha=0.5,label='1x1')
    ax.legend()
    return fig, ax;

if __name__ == "__main__":
    time_a = time.time()
    coords = 'z'
    files = {
        'z1' : f'/rds/general/user/fcm19/home/PhD/photoemission/structures/Cu100_slab_od_tests/Cu100_middle/Cu100_z1/Cu100_od_rerun_{coords}_abs.ome_bin',
        'z3' : f'/rds/general/user/fcm19/home/PhD/photoemission/structures/Cu100_slab_od_tests/Cu100_middle/Cu100_z3/Cu100_od_rerun_{coords}_abs.ome_bin',
        'z5' : f'/rds/general/user/fcm19/home/PhD/photoemission/structures/Cu100_slab_od_tests/Cu100_middle/Cu100_z5/Cu100_od_rerun_{coords}_abs.ome_bin',
    }
    output_file = './Cu100_z_vol_middle.png'
    lower = 1E-15
    upper = 5E2
    fig,ax = make_ome_hist_plot(files=files,coords=coords,lower = lower,upper = upper)
    #ax.set_xlim([2E-1,4E0])
    ax.set_title('Bloch State OMEs for different Cell Volumes')
    plt.tight_layout()
    plt.savefig(output_file,dpi=250)
    time_b = time.time()
    print('It took: ' , (time_b-time_a)/60, 'mins')
