import subprocess
import numpy as np
import os
from pymatgen.core import Lattice
import json
import matplotlib.pyplot as plt
from optados_output_class import *

def read_photonsweep_outputs(seed:str, path = None):
    '''
    This function reads in the output files from a series of optados calculations
    '''
    # create array of the cell volumes, and total energies
    data = {} 

    if path == None:
        path = f'./structures/{seed}/'
    listOfFiles = os.listdir(path)
     # create output classes for each of the output files
    #print(listOfFiles)
    for item in listOfFiles:
        if '.odo' in item:
            #print(item)
            out_obj = OptaDOSOutput(path + item)
            if hasattr(out_obj,'final_state'):
                if out_obj.final_state not in data.keys(): data[out_obj.final_state] = [[out_obj.photon_e,out_obj.qe,out_obj.mte]]
                else: data[out_obj.final_state].append([out_obj.photon_e,out_obj.qe,out_obj.mte])
                data['workfct'] = out_obj.work_fct
    for item in data.keys():
        if item != 'workfct':
            data[item] = np.array(data[item])
            data[item] = data[item][data[item][:,0].argsort()]
    return data;

def make_photonsweep_plots(data:dict,**options):
    plt.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(1,2, figsize = (14,6), dpi = 300)
    if len(data['bloch']) == len(data['free electron']):
        data['mean'] = 0.5*data['bloch']+0.5*data['free electron']
    else: print(f"# of Bloch final state files:{len(data['bloch'])}\n # of Free electron final state files:{len(data['free electron'])}")
    for item in data.keys():
        if item in ['bloch', 'free electron', 'mean']:
            ax[0].plot(data[item][:,0],data[item][:,1], marker = '+', label = f"{item} ({options['temperature']} K)")
            ax[1].plot(data[item][:,0],[x*1000 for x in data[item][:,2]], marker = '+', label =  f"{item} ({options['temperature']} K)")
    
    ax[0].axvline(data['workfct'], ls = '--', c = 'red', label = 'workfunction')
    ax[1].axvline(data['workfct'], ls = '--', c = 'red', label = 'workfunction')
    ax[0].set(xlabel = 'photon energy [eV]', ylabel = 'Quantum efficiency', title = f"{options['title']} quantum efficiency", yscale='log', ylim = [1e-10, 1e-2])
    ax[1].set(xlabel = 'photon energy [eV]', ylabel = 'Mean Transverse Energy (MTE) [meV]', title = f"{options['title']} MTE")
    ax[0].legend()
    ax[1].legend()
    plt.tight_layout()
    return fig,ax;

Cu_surf_data = read_photonsweep_outputs(seed = 'Cu_surf_111', path = '../structures/Cu_surf_111_vic_larger_cell/')
graph_options = {
    'title': 'Cu[111]',
    'temperature': 298,
    'filename' : '../structures/plots/Cu_111_victor_60A_cell_mte_plot_0E.png'
}
#print(Cu_111_data)
fig, ax = make_photonsweep_plots(Cu_surf_data,**graph_options)
plt.savefig(graph_options['filename'],dpi = 250)