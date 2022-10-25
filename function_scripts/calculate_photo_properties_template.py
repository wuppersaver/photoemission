import subprocess
import numpy as np
import os
import sys
from pymatgen.core import Lattice
import json
import matplotlib.pyplot as plt
from optados_output_class import *
import math

def average_potential_from_file(input_file:str, potential = True):
    if potential: factor = 27.211396
    else: factor = 1
    header = []
    with open(input_file, 'r') as f:
       for i in range(12):
            line = next(f)
            header.append(line.strip().split())
            if 'END header:' in line:
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
    # #change file extension
    
    # file_name = input_file.split('/')
    # new_file = file_name[-1]+'.dat'
    # path = ''
    # for item in file_name[:-1]:
    #     path += item + '/'
    # np.savetxt(f'{path}{new_file}',results,delimiter=' ')
    return x_axis*plane_distance, y_axis, cell;

def get_workfct(directory:str=None, file_ending:str = '_optados_photo_sweep.odi', bounds = None,centered:bool = True,mod_odi:bool = True):
    if directory == None:
        directory = f'./structures/' 
    if directory[-1] != '/': directory += '/'
    listOfFiles = os.listdir(directory)
    found = 2
    for item in listOfFiles:
        if '.pot_fmt' in item and not '.dat' in item:
            found = 1
            path = directory + item
            x, potential,cell = average_potential_from_file(path, potential = True)
        if '_fermi.odo' in item:   
            odo_pth = directory + item
            fermi_level = OptaDOSOutput(odo_pth).fermi_e
            found = 0
    if found != 0: 
        file = ['Opatdos File (_fermi.odo)','Potential File (.pot_fmt)']
        raise OSError(2, f'No {file[found-1]} found!')
    indices = [0,0]
    stepsize = x[1] - x[0]
    if not centered and bounds == None: bounds = [int(x[-1]/2)-1,int(x[-1]/2)+1]
    if centered and bounds == None: bounds =[0,5]
    #print(fermi_level.path)
    seed = path.split('/')[-1].split('.')[0]
    for index, item in enumerate(x):
        if abs(item - bounds[0]) <= stepsize: indices[0] = index
        if abs(item - bounds[1]) <= stepsize: indices[1] = index
        if item > (bounds[1]+3*stepsize): break
    vacuum_level = np.mean(potential[indices[0]:indices[1]])
    if mod_odi:             
        for item in listOfFiles:
            if file_ending in item:
                print('Writing work_function=', round(vacuum_level-fermi_level,5),f'eV to {item}')
                subprocess.call(f'sed -i "s/.*work_function.*/work_function : {round(vacuum_level-fermi_level,5)}/" {directory}{item}',shell=True)
    else:
        print('vacuum level = ', round(vacuum_level,5),f' eV, fermi level = ', round(fermi_level,5), ' eV.')
    return;

def get_area_volume(directory:str==None, file_ending:str = '_photo.odi', centered:bool = True,mod_odi:bool = True):
    """This function gets the volume of a slab in a cell by analysing the electron density in the cell. """
    if directory == None:
        directory = f'./structures/' 
    if directory[-1] != '/': directory += '/'
    listOfFiles = os.listdir(directory)
    found = False
    for item in listOfFiles:
        if '.den_fmt' in item and not '.dat' in item:
            path = directory + item
            x, density,cell = average_potential_from_file(path, potential = False)
            found = True
    if not found: raise OSError(2, 'No Density File (.den_fmt) found!')
    #seed = path.split('/')[-1].split('.')[-2]
    area = np.linalg.norm(np.cross(cell.matrix[0],cell.matrix[1]))
    rel_density = density / max(density)
    boundaries = []
    for index, item in enumerate(rel_density):        
        if math.isclose(item,0.01, abs_tol=2e-3):
            boundaries.append(x[index])

    slab_vol =abs(boundaries[-1]-boundaries[0])*area*np.sin(np.deg2rad(cell.alpha))
    if not centered: slab_vol = area*np.linalg.norm(cell.matrix[2]) - slab_vol
    if mod_odi:             
        for item in listOfFiles:
            if file_ending in item:
                print('Writing volume=', round(slab_vol,6), 'A^3 and area=', round(area,5),f'A^2 to {item}')
                subprocess.call(f'sed -i "s/.*surface_area.*/surface_area : {round(area,5)}/" {directory}{item}',shell=True)
                subprocess.call(f'sed -i "s/.*slab_volume.*/slab_volume : {round(slab_vol,6)}/" {directory}{item}',shell=True)
    else:
        print('volume = ', round(slab_vol,6), 'A^3 and area = ', round(area,5),f' A^2')
    return;

if __name__ == "__main__":
    input_path = str(sys.argv[1])
    file_ending = str(sys.argv[2])
    mod_odi = 'yes' in str(sys.argv[3]).lower()
    #input_path = '/rds/general/user/fcm19/home/PhD/photoemission/structures/Cu_surf_100_victor_60A_new/'
    #file_ending = '_optados_photo_sweep.odi',
    get_workfct(directory=input_path, file_ending = file_ending, centered=True,mod_odi=mod_odi)
    get_area_volume(directory=input_path, file_ending = file_ending, centered=True,mod_odi=mod_odi)
