#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 10:50:04 2022

@author: brunocamino
"""
import numpy as np
import subprocess
import ase
import warnings
import json
import os
import math

#from wulffpack import SingleCrystal

from collections import defaultdict

import ase.calculators.castep
import ase.io.castep

from pymatgen.core import Lattice
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import Dos, CompleteDos
from pymatgen.symmetry.kpath import KPathSetyawanCurtarolo
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.eos import EOS

from castep_output_class import *
from optados_output_class import *

import matplotlib.pyplot as plt


def calc_surface_energy(bulk_energy,slab_energy,n_units, area):
    #bulk_energy (float): energy of the bulk
    #slab_energy (float): energy of the slab
    #n_units (int): number of bulk units in the slab
    
    e_surf = 0.5*(slab_energy-(bulk_energy*n_units))/area
    
    return e_surf;

def get_qe(optados_output):
    #This is for you to develop :)
    pass

def generate_input_files(**options):
    task_functions = {
        'Geometry' : generate_castep_input('geometry',**options),
        'Spectral' : generate_castep_input('spectral',**options),
        'OD_Fermi' : generate_optados_input('fermi',**options),
        'Workfunction' : generate_workfct_input('',**options),
        'OD_Photo_Sweep' : generate_optados_input('photo_sweep',**options),
        'BandStructure' : generate_castep_input('bands',**options),
    }
    general = options['general']
    for task in general['tasks']:
        current = task_functions.get(task, lambda: 'Invalid')
    if general['template_scripts'] : generate_scripts(**options)
    return


def generate_scripts(**options):
    blocks = {
        'Geometry' : """            sed -i '0,/.*STATE=.*/s//STATE=geometry_run/' ${calculation[submission]} 
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[geometry]} 
        echo "GeometryOptimization"
        qsub ${calculation[geometry]} 
        exit
    fi
    if [[ $INTERNAL == geometry_run ]]; then    
        sed -i '0,/.*STATE=.*/s//STATE=geometry_cont/' ${calculation[submission]} 
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[geometry]} 
        if ! grep -Fxq "CONTINUATION" ${CASE_IN}_geometry.param; then
            echo "CONTINUATION" | tee -a ${CASE_IN}_geometry.param
        fi
        echo "GeometryOptimization Continued"
        qsub ${calculation[geometry]}
        exit
    fi\n""",
        'Spectral':"""          sed -i '0,/.*STATE=.*/s//STATE=spectral_run/' ${calculation[submission]} 
            sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[spectral]} 
            echo "Spectral"
            qsub ${calculation[spectral]}
            exit
        fi
    fi\n""",
        'OD_Fermi' : """        sed -i '0,/.*STATE=.*/s//STATE=od_fermi_run/' ${calculation[submission]} 
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[optados_all]} 
        echo "OptaDOS misc"
        qsub ${calculation[optados_all]}
        exit
    fi\n""",
        'Workfunction' : """        sed -i '0,/.*STATE=.*/s//STATE=workfct_run/' ${calculation[submission]} 
        sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[work_fct]}
        echo "Setting Workfct and Volume/Area"
        qsub ${calculation[work_fct]}
        exit
    fi\n""",
        'OD_Photo_Sweep': """    sed -i '0,/.*STATE=.*/s//STATE=od_photo_sweep_run/' ${calculation[submission]} 
    sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[optados_photo_sweep]}
    echo "OptaDOS Photoemission Sweep"
    qsub ${calculation[optados_photo_sweep]}
    exit
    fi\n""",
        'BandStructure': """    sed -i '0,/.*STATE=.*/s//STATE=bands_run/' ${calculation[submission]} 
    sed -i '0,/.*CONTINUE=.*/s//CONTINUE=true/'  ${calculation[bands]} 
    echo "Bandstructure"
    qsub ${calculation[bands]}
    exit
    fi\n""",
        'final' : """else
    if [[ $# -gt 0 && $# -lt 2 ]] ; then
        if [ ! -v calculation[$1] ]; then
            echo "Wrong argument was given!!"
            exit 1
        else
            sed -i "0,/.*STATE=.*/s//STATE=${$1}_stop/" ${calculation[submission]} 
            sed -i '0,/.*CONTINUE=.*/s//CONTINUE=false/'  ${calculation[$1]} 
            qsub ${calculation[$1]}
            exit
        fi
    else
        echo "Too many arguments given!!"
        exit 1
    fi
fi\n"""
    }
    task_shorts = {
        'Geometry' : 'geometry',
        'Spectral' : 'spectral',
        'OD_Fermi' : 'od_fermi',
        'Workfunction' : 'workfct',
        'OD_Photo_Sweep' : 'od_photo_sweep',
        'BandStructure' : 'bands',
    }
    tasks = options['general']['tasks']
    directory = options['general']['directory']
    seed = options['general']['seed']
    with open('./templates/template_submission.sh','r') as f:
        lines = f.readlines()
    for index,line in enumerate(lines):
        if 'TEMPLATE' in line:
            lines[index] = line.replace('TEMPLATE',f'{seed}')
    lines.append(blocks[tasks[0]])
    for index,task in enumerate(tasks[1:]):
        lines.append(f"    if [[ $INTERNAL == {task_shorts[tasks[index]]}_success ]]; then\n")
        lines.append(blocks[task])
    lines.append(blocks['final'])
    with open(f'{directory}/{seed}_submission.sh','w') as fw:
        fw.writelines(lines)
    return

def generate_workfct_input(task = None, **options):
    directory = options['general']['directory']
    seed = options['general']['seed']
    lines = []
    with open('./templates/template_workfunction.sh','r') as fr:
        lines = fr.readlines()
    for index,line in enumerate(lines):
        if 'TEMPLATE' in line:
            lines[index] = line.replace('TEMPLATE',f'{seed}')
    with open(f'{directory}/{seed}_workfunction.sh','w') as fw:
        fw.writelines(lines)
    return

def generate_castep_input(task, **options): 
    #based on the CASTEP calculator functionality. For documentation see 
    #https://wiki.fysik.dtu.dk/ase/ase/calculators/castep.html#module-ase.calculators.castep
    
    #calc_struct (ase.Atoms): ASE Atoms object containing the cell and atoms for the calculation
    
    ## Optional Keywords (Must be supplied together)
    ##  directory (str): Directory for generated Files
    ##  label (str): seed or label for CASTEP calculation and files
    ##  task (str): specification of calculation to be performed
    ##  energy_cutoff (int): Energy Cutoff for Plane Waves [eV]
    ##  xc_functional (str): Exchange-Correlation-Functional for the calculation    
    ##  opt_strategy (str): optimization of the calculation (default, speed, 
    ##  memory)
    ##  spin_polarized (str): boolean value, defines if the calculation is spin_polarized, or not
    ##  kpoints (str, list, dict): definition for k-point grid, has dict like functionality
    ##  symmetry (str): boolean value, defines if atoms in cell should be moved to determined symmetry sites
    ##  fix_all_cell (str): boolean value, defines if all the cell parameters should stay fixed during a optimisation run
    ##  continuation (bool): determines if the calculation is a continuation from a previous calc
    tasks = {
        'geometry' : 'GeometryOptimization',
        'spectral' : 'Spectral',
        'single' : 'SinglePoint',
        'bands' : 'BandStructure',
    }
    #print(task)
    castep_task = tasks[task]
    calc_struct = options['structure']
    if not isinstance(calc_struct,ase.atoms.Atoms):
        
        calc_struct = AseAtomsAdaptor().get_atoms(calc_struct)
    general = options['general']
    directory = general['directory']
    seed = general['seed']
    castep = options['castep']
    #print(calc_struct)
    # initialize the calculator instance
    calc = ase.calculators.castep.Castep(check_castep_version = False,keyword_tolerance=3)
    # include interface settings in .param file
    calc._export_settings = False
    
    if options:
        calc._directory = directory
        calc._rename_existing_dir = False
        calc._label = seed
        
        # Define parameter file options
        calc.param.task = castep_task
        if castep_task == 'Spectral': 
            calc.param.spectral_task = castep['spectral_task']
            if castep['calculate_pdos']: calc.param.pdos_calculate_weights = 'TRUE'
        calc.param.xc_functional = castep['xc_functional']
        calc.param.opt_strategy = castep['opt_strategy']
        calc.param.smearing_width = str(castep['smearing_width']) + ' K'
        calc.param.cut_off_energy = str(castep['energy_cutoff']) + ' eV'
        calc.param.elec_energy_tol = str(castep['elec_energy_tol']) + ' eV'
        if castep['extra_bands']: calc.param.perc_extra_bands = '100'
        calc.param.max_scf_cycles = str(castep['max_scf_cycles'])
        if castep['fix_occup']: calc.param.fix_occupancy = 'TRUE'
        if castep['spin_polarized'] : calc.param.spin_polarized = 'TRUE'
        else: calc.param.spin_polarized = 'FALSE'
        if castep['write_potential']: calc.param.write_formatted_potential = 'TRUE'
        if castep['write_density']: calc.param.write_formatted_density = 'TRUE'
        if castep['mixing_scheme'].lower() != 'broydon': calc.param.mixing_scheme = castep['mixing_scheme']
        if castep['continuation']: calc.param.continuation = 'Default'  
        calc.param.num_dump_cycles = 0 # Prevent CASTEP from writing *wvfn* files
        # Define cell file options
        if castep['snap_to_symmetry']: calc.cell.snap_to_symmetry = 'TRUE'
        if castep_task == 'BandStructure':
            band_path = calc_struct.cell.bandpath(castep['bandstruct_path'], density = castep['bandstruct_kpt_dist'])
            print(band_path)
            calc.set_bandpath(bandpath=band_path)
            calc.cell.bs_kpoint_path_spacing = castep['bandstruct_kpt_dist']
        calc.set_kpts(castep['kpoints'])
        if castep_task == 'Spectral': calc.cell.spectral_kpoints_mp_grid = ' '.join([str(x) for x in castep['spectral_kpt_grid']]) 
        if castep_task.lower() == 'geometryoptimization' and castep['fix_all_cell']: calc.cell.fix_all_cell = 'TRUE'
        if castep['generate_symmetry']: calc.cell.symmetry_generate = 'TRUE'
        
    # Prepare atoms and attach them to the current calculator
    calc_struct.calc = calc
   
    #Create input files
    calc_struct.calc.initialize()
    
    # The Cell file has to be modified to have the BS_Kpoint_Path in the right format for CASTEP
    if castep_task == 'BandStructure':
        with open(f"{directory}/{seed}.cell", 'r') as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if 'BS_KPOINT_LIST:' in lines[i]:
                path_list = list(lines[i].replace('BS_KPOINT_LIST: ', '').replace('[','').replace(']', '').replace("'", '').replace('\n','').split(', '))
                path = '%BLOCK BS_KPOINT_PATH\n'
                for point in path_list:
                    path += point + '\n'
                path += '%ENDBLOCK BS_KPOINT_PATH\n'
                lines[i] = path
                break
        with open(f"{directory}/{seed}.cell", 'w') as f:
            for item in lines:
                f.write(item)
    appendices = {
        'geometry' : 'geometry',
        'bands' : 'bands',
        'spectral' : 'spectral',
        'single' : 'single_point',
    }
    if options['general']['template_scripts']:
        lines = []
        with open(f'./templates/template_{appendices[task.lower()]}.sh','r') as fr:
            lines = fr.readlines()
        for index,line in enumerate(lines):
            if 'TEMPLATE' in line:
                lines[index] = line.replace('TEMPLATE',f'{seed}')
        with open(f'{directory}/{seed}_{appendices[task.lower()]}.sh','w') as fw:
            fw.writelines(lines)
    # Rename the Files to allow seamless runs with appropriate shell files
    for file in ['cell','param']: os.replace(f"{directory}/{seed}.{file}",f"{directory}/{seed}_{appendices[task.lower()]}.{file}")
    return;

def generate_optados_input(task,**options):
    directory = options['general']['directory']
    seed = options['general']['seed']
    tasks = {
        'fermi' : 'all',
        'photo' : 'photoemission',
        'photo_sweep' : 'photoemission',
    }
    output = [f"task : {tasks[task]}\n"]
    if tasks[task] != 'photoemission':
        if options['optados']['optados_task'].lower() in ['pdos','all']:
            output.append(f"PDOS : {options['optados']['pdos']}\n")
        output.append(f"broadening : {options['optados']['broadening']}\n")
        output.append(f"iprint : {options['optados']['iprint']}\n")
        output.append(f"efermi : {options['optados']['efermi']}\n")
        output.append(f"dos_spacing : {str(options['optados']['dos_spacing'])}")
    else:
        photo = options['optados']['photo_options']
        for item in photo.keys():
            if item == 'optics_qdir':
                output.append(f"{item} : {photo[item][0]} {photo[item][1]} {photo[item][2]}\n")
                continue
            if item == 'optics_intraband':
                if photo[item]:
                    output.append(f"{item} : TRUE\n")
                continue
            if item == 'hybrid_linear':
                if photo[item]:
                    output.append(f"{item} : TRUE\n")
                continue
            output.append(f"{item} : {photo[item]}\n")
    appendices = {
        'fermi' : 'optados_fermi',
        'photo' : 'optados_photo',
        'photo_sweep' : 'optados_photo_sweep',
        'pdos': 'pdos',
        'dos': 'dos',
    } 
    with open(f"{directory}/{seed}_{appendices[task]}.odi", 'w') as f:
            for line in output:
                f.write(line)
    if options['general']['template_scripts']:
        lines = []
        with open(f'./templates/template_{appendices[task]}.sh','r') as fr:
            lines = fr.readlines()
        if 'sweep' in task:
            sweep_options = options['optados']['sweep_options']
            values = ' '.join(list(sweep_options.values())[1:])
        for index,line in enumerate(lines):
            if 'TEMPLATE' in line:
                lines[index] = line.replace('TEMPLATE',f'{seed}')
            if 'sweep' in task:
                if "sweep_values=seq -f \"%'.5f\"" in line:
                    lines[index] = line.replace('___',values)
                if 's/.*photon_energy.*/photon_energy' in line:
                    lines[index] = line.replace('photon_energy',sweep_options['parameter'])
        with open(f'{directory}/{seed}_{appendices[task]}.sh','w') as fw:
            fw.writelines(lines)
    return;

def get_wulff_fractions(mat_structure:ase.atoms.Atoms, facets_energies : dict):
    oh = SingleCrystal(facets_energies, mat_structure)
    fractions = oh.facet_fractions
    new = defaultdict(list)

    for d in (facets_energies, fractions):
        for key, value in d.items():
            new[key].append(value)
    return new;

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

def create_potential_plot(directory:str=None, bounds = None,centered:bool = True,mod_odi:bool = True):
    if directory == None:
        directory = f'./structures/' 
    if directory[-1] != '/': directory += '/'
    listOfFiles = os.listdir(directory)
    found = [False,False]
    for item in listOfFiles:
        if '.pot_fmt' in item and not '.dat' in item:
            found[0] = True
            path = directory + item
            x, potential,cell = average_potential_from_file(path, potential = True)
        if '_all.odo' in item:   
            odo_pth = directory + item
            fermi_level = OptaDOSOutput(odo_pth).fermi_e
            found[1] = True
        if False not in found: break
    if False in found: 
        file = ['Potential File (.pot_fmt)','Opatdos File (_all.odo)']
        for item in found: 
            if item == False:
                raise OSError(2, f'No {file[found.index[item]]} found!')
    indices = [0,0]
    stepsize = x[1] - x[0]
    if not centered and bounds == None: bounds = [(x[-1]/2)-1,(x[-1]/2)+1]
    if centered and bounds == None: bounds =[0,5]
    #print(fermi_level.path)
    seed = path.split('/')[-1].split('.')[0]
    for index, item in enumerate(x):
        if abs(item - bounds[0]) <= stepsize: indices[0] = index
        if abs(item - bounds[1]) <= stepsize: indices[1] = index
        if item > (bounds[1]+3*stepsize): break
    vacuum_level = np.mean(potential[indices[0]:indices[1]])
    plt.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(1,1, figsize = (12,5), dpi = 300)
    
    ax.hlines(vacuum_level, 0, max(x), ls = '--',colors = 'r', linewidth = 1, label= 'Vacuum Level')
    ax.hlines(fermi_level,0,max(x), ls = '--', colors = 'g', linewidth = 1, label='Fermi Energy')
    ax.vlines(bounds,vacuum_level-1, vacuum_level+1,colors = 'k',linewidth = 1, label = 'bounds')
    ax.plot(x,potential,linewidth = 1, label = 'Potential')
    ax.text(max(x)/2,vacuum_level-2,f'vacuum level = {round(vacuum_level,5)} eV',ha='center')
    ax.text(max(x)/2,fermi_level - 2,f'fermi level = {round(fermi_level,5)} eV',ha='center')
    ax.set_title(f'{seed} - Workfunction W = {round(vacuum_level-fermi_level,5)}')
    ax.set(xlim = (0,max(x)),xlabel = r'Position along c [$\AA$]', ylabel = 'potential [eV]')
    ax.legend(loc='best')
    if mod_odi:             
        for item in listOfFiles:
            if '_photo.odi' in item:
                print('Writing work_function=', round(vacuum_level-fermi_level,5),f'eV to {item}')
                subprocess.call(f'sed -i "s/.*work_function.*/work_function : {round(vacuum_level-fermi_level,5)}/" {directory}{item}',shell=True)

    return fig,ax;

def create_density_plot(directory:str==None, centered:bool = True,mod_odi:bool = True):
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
    seed = path.split('/')[-1].split('.')[-2]
    area = np.linalg.norm(np.cross(cell.matrix[0],cell.matrix[1]))
    rel_density = density / max(density)
    boundaries = []
    for index, item in enumerate(rel_density):        
        if math.isclose(item,0.01, abs_tol=2e-3):
            boundaries.append(x[index])

    slab_vol =abs(boundaries[-1]-boundaries[0])*area*np.sin(np.deg2rad(cell.alpha))
    if not centered: slab_vol = area*np.linalg.norm(cell.matrix[2]) - slab_vol
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
    if mod_odi:             
        for item in listOfFiles:
            if '_photo.odi' in item:
                print('Writing volume=', round(slab_vol,6), 'A^3 and area=', round(area,5),f'A^2 to {item}')
                subprocess.call(f'sed -i "s/.*surface_area.*/surface_area : {round(area,5)}/" {directory}{item}',shell=True)
                subprocess.call(f'sed -i "s/.*slab_volume.*/slab_volume : {round(slab_vol,6)}/" {directory}{item}',shell=True)
    
    return fig,ax;
    
def read_bands2pmg(path:str=None, export = False):
    num_kpoints, num_spin_comps, num_electrons_up, num_electrons_down, num_bands, fermi_energy = 0,0,0,0,0,0
    kpts_coordinates = []
    eigenvalues = {}
    cell = []
    if path == None:
        path = f'./structures/'
    if path[-1] != '/': path += '/'
    listOfFiles = os.listdir(path)
     # create output classes for each of the output files
    for item in listOfFiles:
        if '.bands' in item and '.orig' not in item and '.o2b' not in item:
            seed = item.replace('.bands','')
            with open(f'{path}{item}','r') as f:
                for line in f:
                    line = line.strip()
                    if 'Number of k-points' in line:
                        num_kpoints = int(line.split()[3])
                        num_spin_comps = int(next(f).split()[4])
                        num_electrons = next(f).split()
                        num_electrons_up = float(num_electrons[3])
                        num_bands = int(next(f).split()[3])
                        fermi_energy = float(next(f).split()[5])*27.2113966
                        #print(fermi_energy)

                        kpts_coordinates = np.zeros((num_kpoints,3))
                        eigenvalues[Spin.up] = np.zeros([num_bands, num_kpoints])
                        if num_spin_comps > 1:
                            num_electrons_down = float(num_electrons[4])
                            eigenvalues[Spin.down] = np.zeros([num_bands, num_kpoints])

                        next(f)
                        cell.append([float(x) for x in next(f).split()])
                        cell.append([float(x) for x in next(f).split()])
                        cell.append([float(x) for x in next(f).split()])
                        #print(cell)
                        lattice_obj = Lattice(cell)

                    if line.split()[0] == 'K-point':
                        temp = line.split()
                        index = int(temp[1])-1
                        kpts_coordinates[index] = [float(temp[2]),float(temp[3]),float(temp[4])]
                        next(f)
                        for i in range(num_bands):
                            eigenvalues[Spin.up][i][index] = float(next(f).strip())*27.2113966
                        if num_spin_comps > 1:
                            next(f)
                            for i in range(num_bands):
                                eigenvalues[Spin.down][i][index] = float(next(f).strip())*27.2113966
    for item in listOfFiles: 
        if '_geom.cell' in item or '.cell' in item: 
            cell_item=item
    kpt_path = KPathSetyawanCurtarolo(SpacegroupAnalyzer(read_cell2pmg(f'{path}{cell_item}')).get_primitive_standard_structure()) #Should use the Setyawan-Curtarolo Convention
    high_symm_dict, high_symm_indices = create_label_dict(kpt_path, kpts_coordinates)
    final_kpt_coordinates = np.zeros((num_kpoints+len(high_symm_indices)-2,3))
    final_eigenvalues = {Spin.up : np.zeros([num_bands, num_kpoints+len(high_symm_indices)-2])}
    if num_spin_comps > 1:
        final_eigenvalues = {Spin.down : np.zeros([num_bands, num_kpoints+len(high_symm_indices)-2])}
    
    for i in range(len(high_symm_indices)-1):
        final_kpt_coordinates[high_symm_indices[i]+i:high_symm_indices[i+1]+1+i] = kpts_coordinates[high_symm_indices[i]:high_symm_indices[i+1]+1]
        final_eigenvalues[Spin.up][:,high_symm_indices[i]+i:high_symm_indices[i+1]+1+i] = eigenvalues[Spin.up][:,high_symm_indices[i]:high_symm_indices[i+1]+1]
        if num_spin_comps > 1:
           final_eigenvalues[Spin.down][:,high_symm_indices[i]+i:high_symm_indices[i+1]+1+i] = eigenvalues[Spin.down][:,high_symm_indices[i]:high_symm_indices[i+1]+1]
    new_bandstruct = BandStructureSymmLine(final_kpt_coordinates, final_eigenvalues,lattice_obj.reciprocal_lattice, fermi_energy,high_symm_dict,coords_are_cartesian=False);
    if export:
        with open(f'./structures/band_jsons/{seed}.json', 'w') as f:
            json.dump(new_bandstruct.as_dict(), f)
    return new_bandstruct;

def read_cell2pmg(path:str):
    cell = np.zeros((3,3))
    species = []
    coords = []
    with open(path,'r') as f:
        for line in f:
            line = line.strip()
            if '%BLOCK LATTICE_CART' in line:
                for i in range(3):
                    temp = next(f).strip().split()
                    for j in range(len(temp)):
                        cell[i,j] = float(temp[j])
                lattice_obj = Lattice(cell)
            if '%BLOCK LATTICE_ABC' in line:
                axes = [float(i) for i in next(f).strip().split()]
                angles = [float(i) for i in next(f).strip().split()]
                #print(axes,angles)
                lattice_obj = Lattice.from_parameters(axes[0], axes[1], axes[2], angles[0], angles[1], angles[2])
            if '%BLOCK POSITIONS_ABS' in line:
                while True:
                    temp = next(f).strip().split()
                    if temp[0] == '%ENDBLOCK':
                        break
                    species.append(temp[0])
                    coords.append([float(temp[1]),float(temp[2]),float(temp[3])])
                cartesian = True
                break
            if '%BLOCK POSITIONS_FRAC' in line:
                while True:
                    temp = next(f).strip().split()
                    if temp[0] == '%ENDBLOCK':
                        break
                    species.append(temp[0])
                    coords.append([float(temp[1]),float(temp[2]),float(temp[3])])
                cartesian = False
                break
    return Structure(lattice_obj,species, coords, coords_are_cartesian= cartesian);

def read_geometry_traj_file(path:str = None, seed:str = None):
    ''' This function reads in a castep geometry optimisation trajectory file (.geom) and returns a data dictionary with pymatgen objects, forces and total energies.
    '''
    if path == None:
        if seed == None : raise ValueError('The seed must be specified!')
        path = f'./structures/{seed}/' 
    listOfFiles = os.listdir(path)
    data = {}
    for item in listOfFiles:
        if '.geom' in item:
            with open(path+item) as f:
                for line in f:
                    line = line.strip()
                    if '<-- c' in line:
                        step = int(line.split()[0])
                        line = next(f).strip().split()
                        lattice, data[step] = [], {'total energy' : float(line[0])*27.211386245,'enthalpy' : float(line[1])*27.211386245}
                        for i in range(3): lattice.append([float(x)*0.529177 for x in next(f).strip().split()[0:3]])
                        species, coordinates, line = [], [], next(f).strip()
                        while '<-- R' in line:
                            species.append(line.split()[0])
                            coordinates.append([float(x)*0.529177 for x in line.split()[2:5]])
                            line = next(f).strip()
                        data[step]['structure'] = Structure(lattice, species, coordinates)
                        line, forces = next(f).strip(), []
                        while '<-- F' in line:
                            forces.append([float(x)*51.421 for x in line.split()[2:5]])
                            line = next(f).strip()
                        data[step]['forces'] = forces
    return data;

def create_label_dict(kpath, kpts):
    k_paths = kpath.get_kpoints(line_density=1, coords_are_cartesian = False)
    naming_dict = {}
    labels = {}
    high_symm_indices = []
    for i in range(len(k_paths[0])):
        if k_paths[1][i] != '':
            point = k_paths[1][i]
            temp = [abs(round(k_paths[0][i][0],6)),abs(round(k_paths[0][i][1],6)),abs(round(k_paths[0][i][2],6))]
            naming_dict[point] = tuple(temp)
    for key in naming_dict.keys():
        labels[key] = []
    for i in range(len(kpts)):
        for key, val in naming_dict.items():
            #print(kpts[i],key,val,np.allclose(val,kpts[i]))
            if np.allclose(val,kpts[i]):
                labels[key].append(i)
                high_symm_indices.append(i)
    return naming_dict, high_symm_indices;

def append_symm_line(base_bandstruct, line_to_add):
    base_dict = base_bandstruct.as_dict()
    line_dict = line_to_add.as_dict()

    if base_dict['efermi'] != line_dict['efermi']:
        raise ValueError('Fermi Energies of the two bandstructures are not equal, this will result in faulty concatenation')
    
    base_length = len(base_dict['kpoints'])
    total_kpoints = base_length+len(line_dict['kpoints'])
    num_bands = len(base_dict['bands']['1'])
    
    kpoints = np.zeros((total_kpoints,3))
    eigenvalues = {Spin.up : np.zeros([num_bands, total_kpoints])}
    
    kpoints[:base_length] = base_dict['kpoints']
    kpoints[base_length:] = line_dict['kpoints']
    eigenvalues[Spin.up][:,:base_length] = base_dict['bands']['1']
    eigenvalues[Spin.up][:,base_length:] = line_dict['bands']['1']
    
    if '-1' in base_bandstruct.as_dict()['bands'].keys():
        eigenvalues = {Spin.down : np.zeros([num_bands, total_kpoints])}
        eigenvalues[Spin.down][:,:base_length] = base_dict['bands']['-1']
        eigenvalues[Spin.down][:,base_length:] = line_dict['bands']['-1']
    
    for key in line_dict['labels_dict'].keys():
        base_dict['labels_dict'][key] = base_dict['labels_dict'][key]
    
    return BandStructureSymmLine(kpoints, eigenvalues, base_bandstruct.lattice_rec, base_dict['efermi'], base_dict['labels_dict']);

def plot_dos_optados(path:str=None, plot_up:bool = True, plot_down:bool = False, plot_total:bool = False, xlimit = None, export_json = False):
    energies, up, down, total, max_value = [],[],[],[],[]
    up_total = 0
    down_total = 0
    plt.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(1,1, figsize = (12,5), dpi = 300)
    shift_efermi = True
    if path[-1] != '/':
        path = path + '/'
    if path == None:
        path = f'./structures/' 
    listOfFiles = os.listdir(path)
    for item in listOfFiles:
        if '.adaptive.dat' in item or '.linear.dat' in item:
            seed = item.split('.')[-3]
            with open(f'{path}{item}','r') as f:
                for line in f:
                    line = line.strip()
                    values = line.split()
                    if '#' not in values[0]:
                        energies.append(float(values[0]))
                        up.append(float(values[1]))
                        down.append(float(values[2]))
                        total.append(float(values[1])-float(values[2]))
                        #up_total = float(values[3])
                        #down_total = float(values[4])
        if '.odo' in item:
            with open(f'{path}{item}','r') as o:
                for line in o:
                    line = line.strip()
                    if 'Setting Fermi Energy' in line:
                        temp = next(o).strip().split()
                        efermi = float(temp[6])
                        break
        if '.odi' in item:
            with open(f'{path}{item}','r') as i:
                for line in i:
                    line = line.strip()
                    if 'SET_EFERMI_ZERO' in line:
                        temp = line.split()
                        if temp[2] in ['false','FALSE','False','.false.']:
                            shift_efermi = False
                        break
    
    ax.set(xlim = [energies[0],energies[-1]], xlabel = 'Energy [eV]', ylabel = 'DOS', yticks = [], title = f'DOS Plot for {seed}')
    ax.axhline(0, c = 'gray')
    
    if plot_up:
        max_value.append(max(up))
    if plot_down:
        max_value.append(max(down))
    if plot_total:
        max_value.append(max(total))
    
    if shift_efermi:
        ax.axvline(0, c = '0', ls = '--', alpha = 0.8)
        energies = [item - efermi for item in energies]
        ax.set(xlabel = r'$\mathrm{Energy - E}_{f}$  [eV]', xlim = [energies[0],energies[-1]])
    else:
        ax.axvline(efermi, c = '0', ls = '--', label = 'Fermi Energy', alpha = 0.8)
    
    if plot_up:
        ax.plot(energies, up, label = 'up')
    if plot_down:
        down_flip = [-1*item for item in down]
        ax.plot(energies, down_flip, label = 'down', alpha = 0.5)
    if plot_total:
        ax.plot(energies, total, label = 'total')
    
    ax.set(ylim = [-0.02, max(max_value)+0.1])
    
    if xlimit != None:
        ax.set(xlim = xlimit) 
    
    ax.legend()
    plt.tight_layout()
    total_density = {Spin.up : total}
    if export_json:
        combined_densities =  {Spin.up : up, Spin.down : down}
        if shift_efermi: 
            output_dos = Dos(0, energies,combined_densities)
            total_out = Dos(0, energies, total_density)
        else:
            output_dos = Dos(efermi, energies, combined_densities)
            
        with open(f'./structures/jsons/DOS/{seed}_total.json', 'w') as f:
            json.dump(output_dos.as_dict(), f)
    
    return fig,ax;

def read_proj_dos_optados(path:str=None,  export_json = False,):
    energies, total= [],[]
    columns, projections, column_keys, totals = {}, {}, {}, {Spin.up:[], Spin.down:[]}
    header_string = '#+----------------------------------------------------------------------------+'
    header, values = [],[]
    spin_channels = False
    shifted_efermi = False
    if path == None:
        path = f'./structures/'
    listOfFiles = os.listdir(path)
     # create output classes for each of the output files
    found=[False,False]
    for item in listOfFiles:
        if '_all.odo' in item:
            seed = item.replace('_all.odo','')
            with open(path + item,'r') as g:
                for line in g:
                    if 'Shift energy scale so fermi_energy=0' in line: shifted_efermi = bool(line.split()[7]=='True')
                    if 'Projected Density Of States Calculation' in line:
                        for i in range(6): line = next(g)
                        if 'Fermi energy' in line: efermi = float(line.split()[6])
                        else: print('No fermi energy found, check cursor position!')
            found[0] = True
        if '.pdos.dat' in item:
            with open(path+item,'r') as f:
                for line in f:
                    line_values = line.strip().split()
                    if '#' in line_values[0]:
                        header.append(line_values)
                    else:
                        values.append(line_values)
            found[1] = True
    files=[f'Optados Output',f'pDOS Data File']
    for index,item in enumerate(found):
        if item == False:
            raise  OSError(2, f'No {files[index]} found!')
    for i, item in enumerate(header):
        if 'Column:' in item:
            columns[item[2]] = []
            t = i+2
            while header_string not in header[t]:
                columns[item[2]].append(header[t][1:-1])
                t += 1
    if len(columns['1'][0])>3: spin_channels = True
    #print(columns)
    for key in columns.keys():                                                              # Go through each of the columns read in from the header
        current = columns[key]                                                              # Load the array with the information for the current column, each Atom included has its own array (for example Cu1, Cu2)
        atom = ''
        for item in current:                                                                # Generate string with atoms or sites included in column projector 
            if atom != '':
                atom += ','
            atom += (f'{item[0]}{item[1]}') 
        orbital = current[0][2]                                                             # Read in the orbital projection for this column
        column_keys[key] = [atom,orbital]                                                   # Set the information for that generated column key with information on projector constituents and orbitals
        if atom not in projections.keys():projections[atom] = {}
        if spin_channels:                                                                   # If spin channel is included in the file, initialise initialise projections as {atom: { orbital1: { Spin1: [], ...} , ... } }                    # Initialise atom dict
            if orbital not in projections[atom].keys():projections[atom][orbital] = {}      # Initialise orbital dict
            if spin_channels and current[0][3] in ['Up','Down']: 
                spin = current[0][3]
                if spin not in projections[atom][orbital].keys(): 
                    projections[atom][orbital][spin] = []                                   # initialise spin channel array
                    column_keys[key].append(spin)                                           # Append pin channel to generated column key tuple
        else:                                                                               # If no spin channel is included in the file, initialise projections as {atom: {orbital1: [], ... }}
            if orbital not in projections[atom].keys():projections[atom][orbital] = []
    for item in values:
        item = [float(i) for i in item]
        if not all(abs(elem) == item[1] for elem in item[1:]):
            if not shifted_efermi: energies.append(item[0]-efermi)
            else: energies.append(item[0])
            temp_total, temp_up, temp_down = 0,0,0
            for i in range(len(item[1:])):
                keys = column_keys[str(i+1)]
                if spin_channels:
                    local_dos = float(item[i+1])
                    projections[keys[0]][keys[1]][keys[2]].append(local_dos)
                    if keys[2] == 'Down':
                        temp_total += local_dos * -1
                        temp_down += local_dos
                    else:
                        temp_total += local_dos
                        temp_up += local_dos
                else:
                    projections[keys[0]][keys[1]].append(float(item[i+1]))
                    temp_total = temp_total + float(item[i+1])
            total.append(temp_total)
            totals[Spin.up].append(temp_up)
            totals[Spin.down].append(temp_down)
    if export_json:                                                             # Export currently gives errors due to inconsistencies within pymatgen DOS class                                                             
        cell = read_cell2pmg(path = f"./structures/{seed}/{seed}.cell")
        tot_dos = Dos(energies = energies, densities = totals,efermi = 0)
        proj_dos = CompleteDos(cell, total_dos = tot_dos, pdoss = projections)
        with open(f'./structures/jsons/DOS/{seed}.json', 'w') as f:
            json.dump(proj_dos.as_dict(), f)       
    data = {
        'energies' : energies,
        'total' : totals,
        'projections' : projections,
        'seed' : seed,
        'columns' : columns,
        'spin channels' : spin_channels,

    }
    return data; 

def plot_proj_dos_optados(data:dict=None,plot_up:bool = True, plot_down:bool = False, plot_total:bool = False, xlimit = None):
    plt.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(1,1, figsize = (12,5), dpi = 300)
    if plot_total: ax.plot(data['energies'], data['total'], label = f"{atom}-total", alpha = 0.8,lw = 1)
    for atom in data['projections'].keys():
        for orbital in data['projections'][atom]:
            if data['spin channels']:
                if plot_up: 
                    if not all(abs(elem) == 0.0 for elem in data['projections'][atom][orbital]['Up']): 
                        ax.plot(data['energies'], data['projections'][atom][orbital]['Up'], label = f"{data['columns'][data['columns'].keys()[0]]}-{data['columns'][data['columns'].keys()[-1]]}-{orbital}-Up",alpha = 0.6, lw = 1)
                if plot_down: 
                    if not all(abs(elem) == 0.0 for elem in data['projections'][atom][orbital]['Down']): 
                        ax.plot(data['energies'], data['projections'][atom][orbital]['Down'], label = f"{atom}-{orbital}-Down", alpha = 0.6,lw = 1)
            else:
                keys = data['columns']['1']
                print(keys)
                if not all(abs(elem) == 0.0 for elem in data['projections'][atom][orbital]): 
                    ax.plot(data['energies'], data['projections'][atom][orbital], label = f"{keys[0][0],keys[0][1]}-{keys[-1][0],keys[-1][1]}-{orbital}", alpha = 0.6,lw = 1)
    if xlimit == None: ax.set(xlim = (min(data['energies']),max(data['energies'])))
    else: ax.set(xlim = xlimit)
    plt.legend()
    return fig,ax;

def get_joint_dos_optados(path:str = None, export_json:bool = False):
    if path == None:
        path = f'./structures/'
    if path[-1] != '/': path += '/'
    listOfFiles = os.listdir(path)
     # create output classes for each of the output files
    for item in listOfFiles:
        if '.jadaptive.dat' in item:
            energy,jdos = np.genfromtxt(f'{path}{item}', skip_header=12, dtype =float, unpack=True)
    data = {'energy':energy,'jdos':jdos}
    return data;
    
def plot_joint_dos_optados(data,label, xlimit = [0,12]):
    stepsize = data['energy'][1] - data['energy'][0]
    for index,item in enumerate(data['energy']):
        if xlimit[1] - item < stepsize:
            upper_bound = index
            break
    plt.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(1,1, figsize = (10,5), dpi = 300)
    ax.plot(data['energy'], data['jdos'], label = label,lw = 1)
    ax.set_xlabel('energy [eV]')
    ax.set_ylabel('density of states [electrons/eV]')
    ax.set_ylim([0,data['jdos'][upper_bound]*1.1])
    ax.set_xlim(xlimit)
    return fig,ax;

def add_jdos_line(data,fig,ax,label:str='New jdos'):
    ax.plot(data['energy'],data['jdos'],label = label, lw=1)
    return fig,ax;

def create_slab_layer_convergence(structure, indices, min, max, seed, *cutoff):   
    layers = [i for i in range(min, max+1)]
    castep_opt_template = {
        'directory': f"./structures/{seed}/",
        'label': seed,
        # Param File Instructions
        'task': 'Singlepoint', #Choose: SinglePoint, BandStructure, Spectral 
        'xc_functional': 'PBE',
        'energy_cutoff': 566,
        'opt_strategy': 'Speed',
        'fix_occup' : False,
        'spin_polarized': False,
        'max_scf_cycles': '100',
        'write_potential': False,
        'write_density': False,
        'extra_bands': True,
        'number_extra_bands': 20,
        #Cell File Instructions
        'kpoints': [9,9,9],
    }
    PBS_options = {
        'seed_name': seed,
        'tasks_seeds': [], #Choose one or multiple (carefully!): SinglePoint, BandStructure, Spectral, OptaDOS
        'queue': 'short',
        'castep_version': 'castep_19' #choose from castep_19, castep_18, castep_18_mod
    }
    for i in layers:
        surface = ase.build.surface(lattice = structure, indices =indices, layers = i, vacuum=15, tol=1e-10, periodic=True)
        temp_opt = castep_opt_template
        temp_opt['directory'] = f"./structures/{seed}"
        temp_opt['label'] = f"{seed}_{i}L"
        temp_opt['kpoints'] = get_adjusted_kpts(structure, surface, [9,9,9])
        PBS_options['tasks_seeds'].append(['singlepoint',temp_opt['label']])
        if cutoff != None:
            temp_opt['cutoff'] = cutoff
        generate_castep_input(calc_struct=surface, **temp_opt)
    generate_qsub_file(**PBS_options)
    return;

def create_murnaghan_inputs(seed:str, structure:Structure, cutoff:int, kpoints, dir_path:str = None, min_max:tuple=None):
    # TODO convert this function so that a bash script is created that uses for loops to loop through structures
    '''
    This function takes in a pymatgen structure with a certain cell volume and creates the\\
    input files to run the calculations.
    '''
     # create input path
    if dir_path == None:
        dir_path = f'./structures/{seed}'
    inputs = []
    castep_options = {
        'directory': dir_path,
        'label': seed,
        # Param File Instructions
        'task': 'SinglePoint', #Choose: SinglePoint, BandStructure, Spectral
        'xc_functional': 'PBE',
        'energy_cutoff': cutoff,
        'elec_energy_tol': 1e-8,
        'opt_strategy': 'Speed',
        'fix_occup' : False,
        'mixing_scheme' : 'Pulay', #Choose from Broydon or Pulay
        'smearing_width' : 300, # Smearing width given as a temperature in Kelvin (K)
        'spin_polarized': False,
        'max_scf_cycles': 1000,
        'write_potential': False,
        'write_density': False,
        'extra_bands': True,
        #Cell File Instructions
        'kpoints': kpoints,
        'snap_to_symmetry': False,
        'generate_symmetry': True,
        'fix_all_cell': True,
        'continuation': False,
    }
    PBS_options = {
        'seed_name': seed,
        'tasks_seeds': [['Spectral', seed], ['Optados', seed]], #Choose one or multiple (carefully!): SinglePoint, BandStructure, Spectral, OptaDOS
        'queue': 'short_8_50', # Choose from short_8_50,short_24_100, short_48_100
        'castep_version': 'castep_19' #choose from castep_19, castep_18, castep_18_mod
    }
    # write specific castep options
    castep_options['energy_cutoff'] = cutoff
    castep_options['kpoints'] = kpoints
    castep_options['directory'] = dir_path
    # create the linspace for volumes
    if min_max != None: scaling_factors = np.linspace(min_max[0],min_max[1],min_max[2])
    else: scaling_factors = np.linspace(0.95,1.05,10)
    # scale the given structure lattice and create the input files with the scaled lattices
    initial_volume = structure.lattice.volume
    new_struct = structure
    for factor in scaling_factors:
        factor = round(factor, 3)
        factor_str = str(factor).replace('.','_')
        new_struct.scale_lattice(initial_volume*factor)
        castep_options['label'] = f'{seed}_{factor_str}'
        generate_castep_input(new_struct, **castep_options)
        inputs.append(['SinglePoint',f'{seed}_{factor_str}','castep_18'])
    PBS_options['tasks_seeds'] = inputs
    generate_qsub_file(**PBS_options)
    return;

def read_murnaghan_outputs(seed:str, structure:Structure, path = None):
    '''
    This function reads in the output files from a series of single point calculations\\ 
    and tries to create a pymatgen murnaghan object. Returns this murnaghan object.
    '''
    # create array of the cell volumes, and total energies
    energies, volumes, files = [],[],[]
    if path == None:
        path = f'./structures/{seed}/'
    listOfFiles = os.listdir(path)
     # create output classes for each of the output files
    for item in listOfFiles:
        if '.castep' in item and '.castep_bin' not in item:
            output_temp = CastepOutput(path = path + item)
            files.append([output_temp.seed, output_temp.ks_total_energy, output_temp.structure.lattice.volume])
            energies.append(output_temp.ks_total_energy)
            volumes.append(output_temp.structure.lattice.volume)
    # create pymatgen EOS object with murnaghan
    volumes_sort = np.sort(volumes)
    eos = EOS(eos_name='murnaghan')
    fit = eos.fit(volumes, energies)
    v0 = fit.results['v0']
    min_struct = structure
    min_struct.scale_lattice(v0)
    return fit, min_struct, energies, volumes;

def optados_photon_energy_sweep(seed:str, dir_path:str = None, min_max:tuple=None, **photo_opts):
    '''
    This function takes in a pymatgen structure with a certain cell volume and creates the\\
    input files to run the calculations.
    '''
    OptaDOS_options = {
        'seed_name': seed,
        'optados_task': 'photoemission', # Choose: dos(default), compare_dos, compare_jdos, jdos, pdos, optics, core, all, photoemission
        'broadening': 'adaptive', #Choose: adaptive(default), fixed, linear
        'iprint': '1', #Choose: 1 - bare minimum, 2 - with progress reports, 3 - fulld debug output
        'efermi': 'optados', #Choose: optados - recalculate, file - read from CASTEP file, insulator - count filled bands, float - supplied by user
        'dos_spacing': '0.001', #DOS spacing in unit (default: eV): default - 0.1
        'pdos': 'angular', #Choose: angular, species_ang, species, sites or more detailed descriptions such as: 
        #PDOS : sum:Si1-2(s) - sum of s-chnnls on 2 Si atms (1 proj), 
        #PDOS : Si1;Si2(s) - DOS on Si atom 1 and DOS on s-channel of Si atom 2 (2 proj) 
        'photo_options' : {
            'work_function' : photo_opts['work_function'],
            'surface_area' : photo_opts['surface_area'],
            'slab_volume' : photo_opts['slab_volume'],
            'elec_field' : photo_opts['elec_field'],
            'imfp_const' : 19.0,
            'JDOS_SPACING' : 0.1,
            'JDOS_MAX_ENERGY' : 25,
            'BROADENING' : 'linear',
            'OPTICS_GEOM' : 'unpolar',
            'optics_qdir' : photo_opts['optics_qdir'],
            'photon_energy' : 21.2,
            'linear_smearing' : 0.026,
            'fixed_smearing' :  0.026,
            'optics_intraband' : True,
            'photo_model' : '3step',
            'momentum' : 'crystal',
            'hybrid_linear' : True,
            'temp' : 300,
            'theta_lower' : 59,
            'theta_upper' : 61,
        }
    }
    # create input path
    if dir_path == None:
        dir_path = f'./structures/{seed}'
    # create the linspace for photon energy sweep
    if min_max != None: photon_energies = np.linspace(min_max[0],min_max[1],min_max[2])
    else: photon_energies = np.linspace(OptaDOS_options['photo_options']['work_function']-1.5,OptaDOS_options['photo_options']['work_function']+1.5,10)
    # go through the created photon energies and create the optados inputs
    photon_energies = [round(x,5) for x in photon_energies]
    options = {'optados':OptaDOS_options}
    generate_optados_input(**options)
 
    programs = {'optados': '~/modules_codes/optados/optados.x'}
    output = ["#!/bin/bash  --login\n"]
    output.append(f"#PBS -N {seed}_od_sweep\n")
    output.append('#PBS -l select=1:ncpus=1:mem=50GB\n#PBS -l walltime=02:30:00\n\n')
    output.append("cd $PBS_O_WORKDIR\n\n")
    output.append(f"OPTADOS={programs['optados']}\n\n")
    output.append(f"CASE_IN={seed}\n\n")
    
    energies='energies=(' 
    for energy in photon_energies[:-1]:energies += str(energy) + ' '
    energies += str(photon_energies[-1]) + ')\n\n'
    output.append(energies)
    
    output.append(f"cp ${{CASE_IN}}_photo.odi ${{CASE_IN}}.odi\n\n")

    for state in photo_opts['photo_model']:
        output.append(f"sed -i 's/.*photo_model.*/photo_model : {state}/' ${{CASE_IN}}.odi\n\
CASE_OUT=${{CASE_IN}}_{state}.out\n\n")
        output.append(f"for i in ${{energies[@]}}\n\
do\n\
    sed -i \"s/.*photon_energy.*/photon_energy : $i/\" ${{CASE_IN}}.odi\n\
    $OPTADOS $CASE_IN 2>&1 | tee -a $CASE_OUT \n\
    mv ${{CASE_IN}}.odo ${{CASE_IN}}_${{i}}_{state}.odo\n\
done\n\n")
        
    with open(f"./structures/{seed}/" + f"OptaDOS_sweep_{min_max[0]}-{min_max[1]}.qsub", 'w') as f:
        for line in output:
            f.write(line)
    return;

def read_photonsweep_outputs(path:str = None, seed:str = None):
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
    if options['plottype'] == 'picture': plt.style.use('seaborn-darkgrid')
    if options['plottype'] == 'poster': plt.style.use('seaborn-poster')
    fig, ax = plt.subplots(1,2, figsize = (14,6), dpi = 300)
    if len(data['bloch']) == len(data['free electron']):
        data['mean'] = 0.5*data['bloch']+0.5*data['free electron']
    else: print(f"# of Bloch final state files:{len(data['bloch'])}\n # of Free electron final state files:{len(data['free electron'])}")
    for item in data.keys():
        if item in ['bloch', 'free electron', 'mean']:
            ax[0].plot(data[item][:,0],data[item][:,1], marker = '+', label = f"{item} ({options['temperature']} K)")
            ax[1].plot(data[item][:,0],[x*1000 for x in data[item][:,2]], marker = '+', label =  f"{item} ({options['temperature']} K)")
        if item in ['experimental']:
            ax[1].plot(data[item][:,0],data[item][:,1], c='black' ,marker = '+', label = f"{item} (298 K)")
    ax[0].axvline(data['workfct'], ls = '--', c = 'red', label = 'workfunction')
    ax[1].axvline(data['workfct'], ls = '--', c = 'red', label = 'workfunction')
    ax[0].set(xlabel = 'photon energy [eV]', ylabel = 'Quantum Efficiency (QE)', title = f"{options['title']} Quantum Efficiency", yscale='log', ylim = [1e-10, 1e-2])
    ax[1].set(xlabel = 'photon energy [eV]', ylabel = 'Mean Transverse Energy (MTE) [meV]', title = f"{options['title']} MTE")
    ax[0].legend()
    ax[1].legend()
    plt.tight_layout()
    return fig,ax;

def get_adjusted_kpts(original, new, kpt = [1,1,1]):
    orig_cell  = original.lattice
    new_cell  = new.lattice
    if not all(abs(new_cell.lengths[index] - elem) > 0.001 for index,elem in enumerate(orig_cell.lengths)):
        warnings.warn(f'!WARNING! The cell axes might not be comparable, check the output!')
        print('original\n',orig_cell)
        print('new\n',new_cell)
    new_k = [int(orig_cell.a//new_cell.a*kpt[0]), int(orig_cell.b//new_cell.b*kpt[1]),int(orig_cell.c//new_cell.c*kpt[2])]
    for index, item in enumerate(new_k):
        if item == 0:
            new_k[index] = 1
    return new_k;


def generate_qsub_file(**options):
    # IMPORTANT!!This function is very specific to the machine the generated file is intended for!
    
    # seed_name (str) : seed for input and output files for calculation
    # queue (str) : dictionary key for the submission queue (see queues dict below)
    # program (str) : dict key for the intended calculation program (see programs dict below)
    # bandstructure (bool) : boolean, if bandstructure calculation is performed, it triggers post-processing of *.bands file
    
    
    # !!Programs are specific to Cluster and User!!
    programs = {
        'castep_19':        '~/modules_codes/CASTEP-19.11/obj/linux_x86_64_ifort17--mpi/castep.mpi',
        'castep_18':        '~/modules_codes/CASTEP-18.1/obj/linux_x86_64_ifort17/castep.mpi',
        'castep_18_mod' :   '~/modules_codes/CASTEP-18.1_mod/obj/linux_x86_64_ifort17/castep.mpi',
        'orbitals2bands':   '~/modules_codes/CASTEP-19.11/obj/linux_x86_64_ifort17--serial/orbitals2bands',
        'optados':          '~/modules_codes/optados/optados.x',
        'geom2xyz':         '~/modules_codes/CASTEP-19.11/obj/linux_x86_64_ifort17--serial/geom2xyz'
    }
    pbs = options['pbs']

    output = ["#!/bin/bash  --login\n"]
    output.append(f"#PBS -N {pbs['seed_name']}\n")
    output.append(f"#PBS -l select={pbs['nodes']}:ncpus={pbs['cpus']}:mem={pbs['memory']}GB\n#PBS -l walltime={pbs['walltime']}\n\n")
    output.append("cd $PBS_O_WORKDIR\n\nmodule load mpi intel-suite\n\n")
    
    for task in pbs['tasks_seeds']:
        if task[0].lower()  == 'bandstructure':
            output.append(f"PRGMO2B={programs['orbitals2bands']}\n\n")
        if task[0].lower()  in ['bandstructure','singlepoint','spectral','geometryoptimization']:
            output.append(f"PRGM={programs[task[2]]}\n\n")
            if 'mod' in task[2]: output.append(f"CASE_IN={task[1]}\nCASE_OUT={task[1]}_mod.out\n\n")
            else: output.append(f"CASE_IN={task[1]}\nCASE_OUT={task[1]}.out\n\n")
            output.append("mpiexec $PRGM $CASE_IN 2>&1 | tee $CASE_OUT\n\n")
        if task[0].lower == 'geometryoptimization':
            output.append(f"{programs['geom2xyz']} {pbs['seed_name']}\n")
        if task[0].lower()   == 'bandstructure':
            output.append(f"cp {task[1]}.bands {task[1]}.bands.orig\n\n")
            output.append("$PRGMO2B $CASE_IN 2>&1 | tee -a $CASE_OUT\n\n")
        if task[0].lower()  == 'optados':
            output.append(f"PRGM={programs[task[2]]}\n\n")
            if task[1] != pbs['seed_name']: output.append(f"CASE_IN={pbs['seed_name']}\nCASE_OUT={task[1]}.out\n\n")
            else: output.append(f"CASE_IN={task[1]}\nCASE_OUT={task[1]}_od.out\n\n")
            output.append(f"$PRGM $CASE_IN 2>&1 | tee -a $CASE_OUT\n\n")
        
    with open(f"./structures/{pbs['seed_name']}/{pbs['seed_name']}.qsub", 'w') as f:
        for line in output:
            f.write(line)
    return    
