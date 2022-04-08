#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 10:50:04 2022

@author: brunocamino
"""

from ast import Index
from click import option
import numpy as np
import ase
import re
import warnings
import json

from wulffpack import SingleCrystal

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

def get_energies(seed, path = None):
    scf_lines = []
    if path != None: file = f'{path}/{seed}/{seed}.castep'
    else: file = f'./structures/{seed}/{seed}.castep'
    with open(file,'r') as f:
        for line in f:
            line = line.strip()
            if '*Warning* max. SCF cycles performed' in line:
                warnings.warn(f'!WARNING! The calculation in file {seed}.castep did not converge in the SCF cycle limit.')
            if 'Final energy' in line:
                scf_energy = float(line.split()[4])
            if '-- SCF' in line:
                scf_lines.append(line)
        fermi_energy = float(scf_lines[-2].split()[2])
    return scf_energy,fermi_energy

def generate_qsub_file(**options):
    # IMPORTANT!!This function is very specific to the machine the generated file is intended for!
    
    # seed_name (str) : seed for input and output files for calculation
    # queue (str) : dictionary key for the submission queue (see queues dict below)
    # program (str) : dict key for the intended calculation program (see programs dict below)
    # bandstructure (bool) : boolean, if bandstructure calculation is performed, it triggers post-processing of *.bands file
    
    
    # !!Queue Configurations are specific to Cluster!! SUBJECT TO CHANGE IN APRIL
    queues = {
        'short_8_50': '#PBS -l select=1:ncpus=8:mem=50GB\n#PBS -l walltime=00:30:00\n\n',
        'short_24_100': '#PBS -l select=1:ncpus=24:mem=100GB\n#PBS -l walltime=02:00:00\n\n',
        'short_48_100': '#PBS -l select=1:ncpus=48:mem=100GB\n#PBS -l walltime=02:00:00\n\n',
        
    }
    # !!Programs are specific to Cluster and User!!
    programs = {
        'castep_19':        '/rds/general/user/fcm19/home/modules_codes/CASTEP-19.11/obj/linux_x86_64_ifort17--mpi/castep.mpi',
        'castep_18':        '/rds/general/user/fcm19/home/modules_codes/CASTEP-18.1/obj/linux_x86_64_ifort17/castep.mpi',
        'castep_18_mod' :   '/rds/general/user/fcm19/home/modules_codes/CASTEP-18.1_mod/obj/linux_x86_64_ifort17/castep.mpi',
        'orbitals2bands':   '/rds/general/user/fcm19/home/modules_codes/CASTEP-19.11/obj/linux_x86_64_ifort17--serial/orbitals2bands',
        'optados':          '/rds/general/user/fcm19/home/modules_codes/optados/optados.x'
    }
    
    output = ["#!/bin/bash  --login\n"]
    output.append(f"#PBS -N {options['seed_name']}\n")
    output.append(queues[options['queue']])
    output.append("cd $PBS_O_WORKDIR\n\nmodule load mpi intel-suite\n\n")
    output.append(f"PRGM={programs[options['castep_version']]}\n\n")
    for task in options['tasks_seeds']:
        if task[0].lower()  == 'bandstructure':
            output.append(f"PRGMO2B={programs['orbitals2bands']}\n\n")
        if task[0].lower()  in ['bandstructure','singlepoint','spectral']:
            output.append(f"CASE_IN={task[1]}\nCASE_OUT={task[1]}.out\n\n")
            output.append("mpiexec $PRGM $CASE_IN 2>&1 | tee $CASE_OUT\n\n")
        if task[0].lower()   == 'bandstructure':
            output.append(f"cp {task[1]}.bands {task[1]}.bands.orig\n\n")
            output.append("mpiexec $PRGMO2B $CASE_IN 2>&1 | tee -a $CASE_OUT\n")
        if task[0].lower()  == 'optados':
            output.append(f"OPTADOS={programs['optados']}\n\n")
            output.append(f"CASE_IN={task[1]}\nCASE_OUT={task[1]}.out\n\n")
            output.append(f"mpiexec $OPTADOS $CASE_IN 2>&1 | tee -a $CASE_OUT")
        
    
    with open(f"./structures/{options['seed_name']}/{options['seed_name']}.qsub", 'w') as f:
        for line in output:
            f.write(line)
    return    

def generate_castep_input(calc_struct='hello', **options): 
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
    
    #TODO  
    
    if not isinstance(calc_struct,ase.atoms.Atoms):
        
        calc_struct = AseAtomsAdaptor().get_atoms(calc_struct)

        
    # initialize the calculator instance
    calc = ase.calculators.castep.Castep(check_castep_version = False,keyword_tolerance=3)
    # include interface settings in .param file
    calc._export_settings = False
    
    if options:
        print(options)
        calc._directory = options['directory']
        calc._rename_existing_dir = False
        calc._label = options['label']
        
        # Define parameter file options
        calc.param.task = options['task']
        if options['task'] == 'Spectral': 
            calc.param.spectral_task = options['spectral_task']
            if options['calculate_pdos']: calc.param.pdos_calculate_weights = 'TRUE'
        calc.param.cut_off_energy = str(options['energy_cutoff']) + ' eV'
        calc.param.elec_energy_tol = str(options['elec_energy_tol']) + ' eV'
        calc.param.xc_functional = options['xc_functional']
        calc.param.opt_strategy = options['opt_strategy']
        calc.param.smearing_width = str(options['smearing_width']) + ' K'
        if options['fix_occup']: calc.param.fix_occupancy = 'TRUE'
        if options['spin_polarized'] : calc.param.spin_polarized = 'TRUE'
        else: calc.param.spin_polarized = 'FALSE'
        calc.param.max_scf_cycles = str(options['max_scf_cycles'])
        calc.param.num_dump_cycles = 0 # Prevent CASTEP from writing *wvfn* files
        if options['continuation']: calc.param.continuation = 'Default'  
        if options['write_potential']: calc.param.write_formatted_potential = 'TRUE'
        if options['write_density']: calc.param.write_formatted_density = 'TRUE'
        if options['extra_bands']: calc.param.perc_extra_bands = '100'
        if options['mixing_scheme'].lower() != 'broydon': calc.param.mixing_scheme = options['mixing_scheme']
        # Define cell file options
        if options['snap_to_symmetry']: calc.cell.snap_to_symmetry = 'TRUE'
        if options['task'] == 'BandStructure':
            band_path = calc_struct.cell.bandpath(options['bandstruct_path'], density = options['bandstruct_kpt_dist'])
            calc.set_bandpath(bandpath=band_path)
            calc.cell.bs_kpoint_path_spacing = options['bandstruct_kpt_dist']
        calc.set_kpts(options['kpoints'])
        if options['task'] == 'Spectral': calc.cell.spectral_kpoints_mp_grid = f"{options['spectral_kpt_grid'][0]} {options['spectral_kpt_grid'][1]} {options['spectral_kpt_grid'][2]}"
        if options['fix_all_cell']: calc.cell.fix_all_cell = 'TRUE'
        if options['generate_symmetry']: calc.cell.symmetry_generate = 'TRUE'
        
    # Prepare atoms and attach them to the current calculator
    calc_struct.calc = calc
   
    #Create input files
    calc_struct.calc.initialize()
    
    # The Cell file has to be modified to have the BS_Kpoint_Path in the right format for CASTEP
    if options['task'] == 'BandStructure':
        with open(f"./structures/{options['label']}/{options['label']}.cell", 'r') as f:
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
        with open(f"./structures/{options['label']}/{options['label']}.cell", 'w') as f:
            for item in lines:
                f.write(item)

    return;

def generate_optados_input(options, **photo):
    output = [f"task : {options['optados_task']}\n"]
    if options['optados_task'].lower() != 'photoemission':
        if options['optados_task'].lower() == 'pdos':
            output.append(f"pdos : {options['pdos']}\n")
        output.append(f"broadening : {options['broadening']}\n")
        output.append(f"iprint : {options['iprint']}\n")
        output.append(f"efermi : {options['efermi']}\n")
        output.append(f"dos_spacing : {str(options['dos_spacing'])}")
    else:
        for item in photo.keys():
            if item == 'optics_qdir':
                output.append(f"{item} : {photo[item][0]} {photo[item][1]}{photo[item][2]}")
                pass
            if item == 'optics_intraband':
                if photo[item]:
                    output.append(f"{item} : TRUE")
                pass
            if item == 'hybrid_linear':
                if photo[item]:
                    output.append(f"{item} : TRUE")
                pass
            output.append(f"{item} : {photo[item]}")
    with open(f"./structures/{options['seed_name']}/{options['seed_name']}.odi", 'w') as f:
        for line in output:
            f.write(line)
    return;

def get_wulff_fractions(mat_structure:ase.atoms.Atoms, facets_energies : dict):
    oh = SingleCrystal(facets_energies, mat_structure)
    fractions = oh.facet_fractions
    new = defaultdict(list)

    for d in (facets_energies, fractions):
        for key, value in d.items():
            new[key].append(value)
    return new;

def average_potential_from_file(input_file:str):

    spin = np.genfromtxt(input_file, skip_header=7, delimiter=None, invalid_raise=False)
    
    spin_value = int(spin[0])
    
    x_y_z = np.genfromtxt(input_file, delimiter=None, skip_header=8, invalid_raise=False)
    #extracts x, y, z axes from slab and creates array
    
    number_of_points_per_plane = int(x_y_z[0]*x_y_z[1])*spin_value
    #gets number of points per plane in z-axis of slab
    
    number_of_planes = int(x_y_z[2]+1)
    #needs '+1' since when used later to generate x-axis of graph np.arange() is exclusive of end of range
    
    data_array = np.genfromtxt(input_file, skip_header=11, delimiter=None)
    
    sorted_data = data_array[data_array[:,2].argsort()] #sorts data by third column
    
    energy_column = sorted_data[:,3] #extracts energy column 
    
    total_energy_of_each_plane = np.add.reduceat(energy_column, range(0,len(energy_column), number_of_points_per_plane))
    #adds up energies in chunks corresponding to how many planes there are
    
    mean_energy_of_plane = total_energy_of_each_plane / number_of_points_per_plane
    
    x_axis = np.arange(1, number_of_planes)
    y_axis = mean_energy_of_plane * 27.2114
    
    results = np.c_[x_axis,y_axis]
    #change file extension
    
    file_name = input_file.split('.')
    new_file = file_name[0]+'.dat'
    
    np.savetxt(new_file,results,delimiter=' ')
    return;
    
def read_bands2pmg(seed:str, export = False):
    num_kpoints, num_spin_comps, num_electrons_up, num_electrons_down, num_bands, fermi_energy = 0,0,0,0,0,0
    kpts_coordinates = []
    eigenvalues = {}
    cell = []
    with open(f'./structures/{seed}/{seed}.bands','r') as f:
        for line in f:
            line = line.strip()
            if 'Number of k-points' in line:
                num_kpoints = int(line.split()[3])
                num_spin_comps = int(next(f).split()[4])
                num_electrons = next(f).split()
                num_electrons_up = float(num_electrons[3])
                num_bands = int(next(f).split()[3])
                fermi_energy = float(next(f).split()[5])*27.2113966

                kpts_coordinates = np.zeros((num_kpoints,3))
                eigenvalues[Spin.up] = np.zeros([num_bands, num_kpoints])
                if num_spin_comps > 1:
                    num_electrons_down = float(num_electrons[4])
                    eigenvalues[Spin.down] = np.zeros([num_bands, num_kpoints])
                
                next(f)
                cell.append([float(x) for x in next(f).split()])
                cell.append([float(x) for x in next(f).split()])
                cell.append([float(x) for x in next(f).split()])
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
    kpt_path = KPathSetyawanCurtarolo(SpacegroupAnalyzer(read_cell2pmg(f'./structures/{seed}/{seed}.cell')).get_primitive_standard_structure()) #Should use the Setyawan-Curtarolo Convention
    
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
    return new_bandstruct

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
                lattice_obj = Lattice.from_parameters(axes[0], axes[1], axes[2], angles[0], angles[1], angles[1])
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

def create_label_dict(kpath, kpts):
    k_paths = kpath.get_kpoints(line_density=1, coords_are_cartesian = False)
    naming_dict = {}
    labels = {}
    high_symm_indices = []
    for i in range(len(k_paths[0])):
        if k_paths[1][i] != '':
            point = k_paths[1][i]
            temp = [round(k_paths[0][i][0],3),round(k_paths[0][i][1],3),round(k_paths[0][i][2],3)]
            for i in range(3):
                if temp[i] == -0.0:
                    temp[i] = 0.0
            naming_dict[point] = tuple(temp)
    for key in naming_dict.keys():
        labels[key] = []
    for i in range(len(kpts)):
        for key, val in naming_dict.items():
            if np.array_equal(val,kpts[i]):
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

def plot_dos_optados(seed:str, plot_up:bool = True, plot_down:bool = False, plot_total:bool = False, xlimit = None, export_json = False):
    energies, up, down, total, max_value = [],[],[],[],[]
    up_total = 0
    down_total = 0
    plt.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(1,1, figsize = (12,5), dpi = 300)
    shift_efermi = True
    with open(f'./structures/{seed}/{seed}.adaptive.dat','r') as f:
        for line in f:
            line = line.strip()
            values = line.split()
            if '#' not in values[0]:
                energies.append(float(values[0]))
                up.append(float(values[1]))
                down.append(float(values[2]))
                total.append(float(values[1])-float(values[2]))
                up_total = float(values[3])
                down_total = float(values[4])
    with open(f"./structures/{seed}/{seed}.odo", 'r') as o:
        for line in o:
            line = line.strip()
            if 'Setting Fermi Energy' in line:
                temp = next(o).strip().split()
                efermi = float(temp[6])
                break
    with open(f"./structures/{seed}/{seed}.odi", "r") as i:
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

def plot_proj_dos_optados(seed:str, plot_up:bool = True, plot_down:bool = False, plot_total:bool = False, xlimit = None, export_json = False):
    energies, total= [],[]
    columns, projections, column_keys, totals = {}, {}, {}, {Spin.up:[], Spin.down:[]}
    plt.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(1,1, figsize = (12,5), dpi = 300)
    header_string = '#+----------------------------------------------------------------------------+'
    header, values = [],[]
    spin_channels = False
    with open(f'./structures/{seed}/{seed}.pdos.dat','r') as f:
        for line in f:
            line_values = line.strip().split()
            if '#' in line_values[0]:
                header.append(line_values)
            else:
                values.append(line_values)
    for i, item in enumerate(header):
        if 'Column:' in item:
            columns[item[2]] = []
            t = i+2
            while header_string not in header[t]:
                columns[item[2]].append(header[t][1:-1])
                t += 1
    if len(columns['1'][0])>3: spin_channels = True
    print(columns)
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
            energies.append(item[0])
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
    if plot_total: ax.plot(energies, total, label = f"{atom}-total", alpha = 0.8,lw = 1)
    for atom in projections.keys():
        for orbital in projections[atom]:
            if spin_channels:
                if plot_up: 
                    if not all(abs(elem) == 0.0 for elem in projections[atom][orbital]['Up']): 
                        ax.plot(energies, projections[atom][orbital]['Up'], label = f"{atom}-{orbital}-Up",alpha = 0.6, lw = 1)
                if plot_down: 
                    if not all(abs(elem) == 0.0 for elem in projections[atom][orbital]['Down']): 
                        ax.plot(energies, projections[atom][orbital]['Down'], label = f"{atom}-{orbital}-Down", alpha = 0.6,lw = 1)
            else:
                if not all(abs(elem) == 0.0 for elem in projections[atom][orbital]): ax.plot(energies, projections[atom][orbital], label = f"{atom}-{orbital}", alpha = 0.6,lw = 1)
    if xlimit == None: ax.set(xlim = (min(energies),max(energies)))
    else: ax.set(xlim = xlimit)
    plt.legend()
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

