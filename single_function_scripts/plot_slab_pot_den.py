import subprocess
import numpy as np
import os
from pymatgen.core import Lattice
import json
import matplotlib.pyplot as plt
from optados_output_class import *

def average_potential_from_file(input_file:str, potential = True):
    if potential: factor = 27.211396
    else: factor = 1
    header = []
    with open(input_file, 'r') as f:
        while True:
            line = next(f)
            header.append(line.strip().split())
            if 'END' in line:
                break

    spin_value = int(header[7][0])
    
    lengths = [float(x[5]) for x in header[3:6]]
    lattice = list(map(lambda sub: list(map(float, sub)), [x[0:3] for x in header[3:6]])) 
    cell = Lattice(lattice)
    
    number_of_points_per_plane = int(header[8][0])*int(header[8][1])*spin_value
    #gets number of points per plane in z-axis of slab
    number_of_planes = int(header[8][2])+1
    plane_distance = lengths[2]/number_of_planes
    
    #needs '+1' since when used later to generate x-axis of graph np.arange() is exclusive of end of range
    
    data_array = np.genfromtxt(input_file, skip_header=11, delimiter=None)
    
    sorted_data = data_array[data_array[:,2].argsort()] #sorts data by third column
    
    energy_column = sorted_data[:,3] #extracts energy column 
    
    total_energy_of_each_plane = np.add.reduceat(energy_column, range(0,len(energy_column), number_of_points_per_plane))
    #adds up energies in chunks corresponding to how many planes there are
    
    mean_energy_of_plane = total_energy_of_each_plane / number_of_points_per_plane
    
    x_axis = np.arange(1, number_of_planes, dtype=np.int16)
    y_axis = mean_energy_of_plane * factor
    
    results = np.c_[x_axis,y_axis]
    #change file extension
    
    file_name = input_file.split('/')
    new_file = file_name[-1]+'.dat'
    path = ''
    for item in file_name[:-1]:
        path += item + '/'
    np.savetxt(f'{path}{new_file}',results,delimiter=' ')
    return x_axis*plane_distance, y_axis, cell;

def create_potential_plot(file_path:str, bounds = [0,5.4]):
    x, potential,cell = average_potential_from_file(file_path, potential = True)
    indices = [0,0]
    stepsize = x[1] - x[0]
    fermi_level = OptaDOSOutput(file_path.replace('.pot_fmt','_fermi.odo')).fermi_e
    #print(fermi_level.path)
    seed = file_path.split('/')[-1].split('.')[0]
    for index, item in enumerate(x):
        if abs(item - bounds[0]) <= stepsize: indices[0] = index
        if abs(item - bounds[1]) <= stepsize: indices[1] = index
        if item > (bounds[1]+3*stepsize): break
    mean = np.mean(potential[indices[0]:indices[1]])
    plt.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(1,1, figsize = (12,5), dpi = 300)
    
    ax.hlines(mean, 0, max(x), ls = '--',colors = 'r', linewidth = 1, label= 'Vacuum Level')
    ax.hlines(fermi_level,0,max(x), ls = '--', colors = 'g', linewidth = 1, label='Fermi Energy')
    ax.plot(x,potential,linewidth = 1, label = 'Potential')
    ax.text(max(x)/2,mean-2,f'vacuum level = {round(mean,5)} eV',ha='center')
    ax.text(max(x)/2,fermi_level - 2,f'fermi level = {round(fermi_level,5)} eV',ha='center')
    ax.set_title(f'{seed} - Workfunction W = {round(mean-fermi_level,5)}')
    ax.set(xlim = (0,max(x)),xlabel = r'Position along c [$\AA$]', ylabel = 'potential [eV]')
    ax.legend(loc='best')

    return fig,ax;

def create_density_plot(file_path:str):
    x, density,cell = average_potential_from_file(file_path, potential = False)
    seed = file_path.split('/')[-1].split('.')[-2]
    area = np.linalg.norm(np.cross(cell.matrix[0],cell.matrix[1]))
    rel_density = density / max(density)
    boundaries = []
    for index, item in enumerate(rel_density):
        if abs(item-0.01) < 0.002:
            boundaries.append(x[index])

    slab_vol =abs(boundaries[-1]-boundaries[0])*area*np.sin(np.deg2rad(cell.alpha))

    plt.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(1,1, figsize = (12,5), dpi = 300)

    ax.plot(x,rel_density,marker = '+', label = 'charge density')
    ax.vlines(boundaries,-0.1,1,colors='r', ls = '--', linewidth = 1, label= 'slab boundaries')
    #ax.text(max(x) - 12, 0.5, fr'Area = {round(area,5)} $\AA^2$')
    #ax.text(max(x) - 12, 0.4, fr'Slab Volume = {round(slab_vol,5)} $\AA^3$')
    ax.set_title(fr'{seed} Area = {round(area,5)} $\AA^2$ - Slab Volume = {round(slab_vol,5)} $\AA^3$')
    ax.legend(loc = 'best')
    #ax.set(xlim = (10,15),ylim = (-0.001,0.015))# max(rel_density)+0.05))
    ax.set(xlim = (0,max(x)),ylim = (-0.1,max(rel_density)+0.05),xlabel = r'Position along c [$\AA$]', ylabel = 'scaled electronic density')
    plt.tight_layout()
    print(round(slab_vol,5))
    return fig,ax;

seed = 'Cu_surf_111'

fig,ax = create_density_plot(file_path=f'../structures/{seed}/{seed}.den_fmt')
plt.tight_layout()
plt.savefig(f'../structures/plots/{seed}_victor_large_16L_potential_plot.png')

fig,ax = create_potential_plot(file_path=f'../structures/{seed}/{seed}.pot_fmt', bounds= [25,27])
plt.tight_layout()
plt.savefig(f'../structures/plots/{seed}_victor_large_cell_16L_potential_plot.png')