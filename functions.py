#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 10:50:04 2022

@author: brunocamino
"""

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
from pymatgen.symmetry.kpath import KPathSetyawanCurtarolo
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


import matplotlib.lines as mlines
import math
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

def get_energies(seed, different_path = False, path = ''):
    lines = []
    scf_lines = []
    scf_energy = 0  
    if different_path:
        with open(f'./{path}{seed}/{seed}.castep','r') as f:
                 lines = f.readlines()
    else:
        with open(f'./structures/{seed}/{seed}.castep','r') as f:
                 lines = f.readlines()
    for line in lines:
        if '*Warning* max. SCF cycles performed' in line:
            warnings.warn(f'!WARNING! The calculation in file {seed}.castep did not converge in the SCF cycle limit.')
        if 'Final energy' in line:
            scf_energy = float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)",line)[0])
        if '-- SCF' in line:
            scf_lines.append(line)
    fermi_energy = float(re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?',scf_lines[-2])[2])
    return scf_energy,fermi_energy

def generate_qsub_file(**options):
    # IMPORTANT!!This function is very specific to the machine the generated file is intended for!
    
    # seed_name (str) : seed for input and output files for calculation
    # queue (str) : dictionary key for the submission queue (see queues dict below)
    # program (str) : dict key for the intended calculation program (see programs dict below)
    # bandstructure (bool) : boolean, if bandstructure calculation is performed, it triggers post-processing of *.bands file
    
    
    # !!Queue Configurations are specific to Cluster!! SUBJECT TO CHANGE IN APRIL
    queues = {
        'short': '#PBS -l select=1:ncpus=48:mem=100GB\n#PBS -l walltime=02:00:00\n\n',
        'debug': '#PBS -l select=1:ncpus=8:mem=90GB\n#PBS -l walltime=00:30:00\n\n'
    }
    # !!Programs are specific to Cluster and User!!
    programs = {
        'castep_19': '/rds/general/user/fcm19/home/modules_codes/CASTEP-19.11/obj/linux_x86_64_ifort17--mpi/castep.mpi',
        'castep_18': '/rds/general/user/fcm19/home/modules_codes/CASTEP-18.1/obj/linux_x86_64_ifort17/castep.mpi',
        'castep_18_mod' : '/rds/general/user/fcm19/home/modules_codes/CASTEP-18.1_mod/obj/linux_x86_64_ifort17/castep.mpi',
        'orbitals2bands': '/rds/general/user/fcm19/home/modules_codes/CASTEP-19.11/obj/linux_x86_64_ifort17--serial/orbitals2bands',
        'optados': '/rds/general/user/fcm19/home/modules_codes/optados/optados.x'
    }
    
    output = ["#!/bin/bash  --login\n"]
    output.append(f"#PBS -N {options['seed_name']}\n")
    output.append(queues[options['queue']])
    output.append("cd $PBS_O_WORKDIR\n\nmodule load mpi intel-suite\n\n")
    output.append(f"PRGM={programs[options['castep_version']]}\n\n")
    for task in options['tasks']:
        if task == 'Bandstructure':
            output.append(f"PRGMO2B={programs['orbitals2bands']}\n\n")
        if task in ['Bandstructure','SinglePoint','Spectral']:
            output.append(f"CASE_IN={options['seed_name']}\nCASE_OUT={options['seed_name']}.out\n\n")
            output.append("mpiexec $PRGM $CASE_IN 2>&1 | tee $CASE_OUT\n\n")
        if task == 'Bandstructure':
            output.append(f"cp {options['seed_name']}.bands {options['seed_name']}.bands.orig\n\n")
            output.append("mpiexec $PRGMO2B $CASE_IN 2>&1 | tee -a $CASE_OUT\n")
        if task == 'Optados':
            output.append(f"OPTADOS={programs['optados']}\n\n")
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
        raise TypeError('An ASE Atoms object has to be given!')        
        
    # initialize the calculator instance
    calc = ase.calculators.castep.Castep(check_castep_version = False,keyword_tolerance=3)
    # include interface settings in .param file
    calc._export_settings = False

    if options:
        calc._directory = options['directory']
        calc._rename_existing_dir = False
        calc._label = options['label']
        
        # Define parameter file options
        calc.param.task = options['task']
        if options['task'] == 'Spectral': 
            calc.param.spectral_task = options['spectral_task']
            if options['calculate_pdos']:
                calc.param.pdos_calculate_weights = 'TRUE'
        calc.param.cut_off_energy = str(options['energy_cutoff']) + ' eV'
        calc.param.xc_functional = options['xc_functional']
        calc.param.opt_strategy = options['opt_strategy']
        if options['fix_occup']: calc.param.fix_occupancy = 'TRUE'
        if options['spin_polarized'] : calc.param.spin_polarized = 'TRUE'
        calc.param.max_scf_cycles = options['max_scf_cycles']
        calc.param.num_dump_cycles = 0 # Prevent CASTEP from writing *wvfn* files
        if options['continuation']: calc.param.continuation = 'Default'  
        if options['write_potential']: calc.param.write_formatted_potential = 'TRUE'
        if options['write_density']: calc.param.write_formatted_density = 'TRUE'
        if options['extra_bands']: calc.param.nextra_bands = options['number_extra_bands']
        
        # Define cell file options
        if options['snap_to_symmetry']: calc.cell.snap_to_symmetry = 'TRUE'
        if options['task'] == 'BandStructure':
            band_path = calc_struct.cell.bandpath(options['bandstruct_path'], density = options['bandstruct_kpt_dist'])
            calc.set_bandpath(bandpath=band_path)
            calc.cell.bs_kpoint_path_spacing = options['bandstruct_kpt_dist']
        calc.set_kpts(options['kpoints'])
        if options['task'] == 'Spectral': calc.cell.spectral_kpoints_mp_grid = f"{options['spectral_kpt_grid'][0]} {options['spectral_kpt_grid'][1]} {options['spectral_kpt_grid'][2]}"
        if options['fix_all_cell']: calc.cell.fix_all_cell = 'TRUE'
       
        
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

def generate_optados_input(options):
    output = [f"task : {options['optados_task']}\n"]
    output.append(f"broadening : {options['broadening']}\n")
    output.append(f"iprint : {options['iprint']}\n")
    output.append(f"efermi : {options['efermi']}\n")
    output.append(f"dos_spacing : {str(options['dos_spacing'])}")
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
                
            if '%BLOCK POSITIONS_ABS' in line:
                while True:
                    temp = next(f).strip().split()
                    if temp[0] == '%ENDBLOCK':
                        break
                    species.append(temp[0])
                    coords.append([float(temp[1]),float(temp[2]),float(temp[3])])
                break
                
    lattice_obj = Lattice(cell)
    return Structure(lattice_obj,species, coords);

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

def plot_dos_optados(seed:str, plot_up:bool = True, plot_down:bool = False, plot_total:bool = False):
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
                if temp[2] == 'false':
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
    ax.legend()
    plt.tight_layout()
    return fig,ax;

#Everything below this is to plot bandstructure plots with functions adapted from the pymatgen module

def get_bs_plot(
        bandstructure,
        zero_to_efermi=True,
        ylim=None,
        smooth=False,
        vbm_cbm_marker=False,
        smooth_tol=0,
        smooth_k=3,
        smooth_np=100,
        bs_labels=[],
    ):
        plt.style.use('seaborn-darkgrid')
        plt.rcParams["font.family"] = "serif"
        pretty = pretty_plot(12, 8)
        bs_array = [bandstructure]

        if isinstance(smooth, bool):
            smooth = [smooth] * len(bs_array)

        handles = []
        vbm_min, cbm_max = [], []

        colors = list(pretty.rcParams["axes.prop_cycle"].by_key().values())[0]
        for ibs, bs in enumerate(bs_array):

            # set first bs in the list as ref for rescaling the distances of the other bands
            bs_ref = bs_array[0] if len(bs_array) > 1 and ibs > 0 else None

            if smooth[ibs]:
                # interpolation works good on short segments like branches
                data = bs_plot_data(zero_to_efermi, bs, bs_ref, split_branches=True)
            else:
                data = bs_plot_data(zero_to_efermi, bs, bs_ref, split_branches=False)

            # remember if one bs is a metal for setting the ylim later
            one_is_metal = False
            if not one_is_metal and data["is_metal"]:
                one_is_metal = data["is_metal"]

            # remember all the cbm and vbm for setting the ylim later
            if not data["is_metal"]:
                cbm_max.append(data["cbm"][0][1])
                vbm_min.append(data["vbm"][0][1])

            for sp in bs.bands.keys():
                ls = "-" if str(sp) == "1" else "--"

                if bs_labels != []:
                    bs_label = f"{bs_labels[ibs]} {sp.name}"
                else:
                    bs_label = f"Band {ibs} {sp.name}"

                handles.append(mlines.Line2D([], [], lw=2, ls=ls, color=colors[ibs], label=bs_label))

                distances, energies = data["distances"], data["energy"][str(sp)]

                for dist, ene in zip(distances, energies):
                    pretty.plot(dist, ene.T, ls=ls)

            # plot markers for vbm and cbm
            if vbm_cbm_marker:
                for cbm in data["cbm"]:
                    pretty.scatter(cbm[0], cbm[1], color="r", marker="o", s=100)
                for vbm in data["vbm"]:
                    pretty.scatter(vbm[0], vbm[1], color="g", marker="o", s=100)

            # Draw Fermi energy, only if not the zero
            if not zero_to_efermi:
                ef = bs.efermi
                pretty.axhline(ef, lw=2, ls="-.", color=colors[ibs])

        # defaults for ylim
        e_min = -4
        e_max = 4
        if one_is_metal:
            e_min = -10
            e_max = 10

        if ylim is None:
            if zero_to_efermi:
                if one_is_metal:
                    # Plot A Metal
                    pretty.ylim(e_min, e_max)
                else:
                    pretty.ylim(e_min, max(cbm_max) + e_max)
            else:
                all_efermi = [b.efermi for b in bs_array]
                ll = min([min(vbm_min), min(all_efermi)])
                hh = max([max(cbm_max), max(all_efermi)])
                pretty.ylim(ll + e_min, hh + e_max)
        else:
            pretty.ylim(ylim)
        pretty = maketicks(bs, pretty)

        # Main X and Y Labels
        #plt.xlabel(r"$\mathrm{Wave\ Vector}$", fontsize=30)
        pretty.xlabel('')
        ylabel = r"$\mathrm{E\ -\ E_f\ (eV)}$" if zero_to_efermi else r"$\mathrm{Energy\ (eV)}$"
        pretty.ylabel(ylabel, fontsize=20)

        # X range (K)
        # last distance point
        x_max = data["distances"][-1][-1]
        pretty.xlim(0, x_max)

        #plt.legend(handles=handles)

        pretty.tight_layout()

        # auto tight_layout when resizing or pressing t
        def fix_layout(event):
            if (event.name == "key_press_event" and event.key == "t") or event.name == "resize_event":
                pretty.gcf().tight_layout()
                pretty.gcf().canvas.draw()

        pretty.gcf().canvas.mpl_connect("key_press_event", fix_layout)
        pretty.gcf().canvas.mpl_connect("resize_event", fix_layout)

        return pretty

def pretty_plot(width=8, height=None, plt=None, dpi=300, color_cycle=("qualitative", "Set1_9")):
    """
    Provides a publication quality plot, with nice defaults for font sizes etc.

    Args:
        width (float): Width of plot in inches. Defaults to 8in.
        height (float): Height of plot in inches. Defaults to width * golden
            ratio.
        plt (matplotlib.pyplot): If plt is supplied, changes will be made to an
            existing plot. Otherwise, a new plot will be created.
        dpi (int): Sets dot per inch for figure. Defaults to 300.
        color_cycle (tuple): Set the color cycle for new plots to one of the
            color sets in palettable. Defaults to a qualitative Set1_9.

    Returns:
        Matplotlib plot object with properly sized fonts.
    """
    ticksize = int(width * 2.5)

    golden_ratio = (math.sqrt(5) - 1) / 2

    if not height:
        height = int(width * golden_ratio)

    if plt is None:
        import importlib

        import matplotlib.pyplot as plt

        mod = importlib.import_module(f"palettable.colorbrewer.{color_cycle[0]}")
        colors = getattr(mod, color_cycle[1]).mpl_colors
        from cycler import cycler

        plt.figure(figsize=(width, height), facecolor="w", dpi=dpi)
        ax = plt.gca()
        ax.set_prop_cycle(cycler("color", colors))
    else:
        fig = plt.gcf()
        fig.set_size_inches(width, height)
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)

    ax = plt.gca()
    ax.set_title(ax.get_title(), size=width * 4)

    labelsize = int(width * 3)

    ax.set_xlabel(ax.get_xlabel(), size=labelsize)
    ax.set_ylabel(ax.get_ylabel(), size=labelsize)

    return plt

def bs_plot_data(zero_to_efermi=True,bs = None, bs_ref=None, split_branches=True):
        """
        Get the data nicely formatted for a plot

        Args:
            zero_to_efermi: Automatically subtract off the Fermi energy from the
                eigenvalues and plot.
            bs: the bandstructure to get the data from. If not provided, the first
                one in the self._bs list will be used.
            bs_ref: is the bandstructure of reference when a rescale of the distances
                is need to plot multiple bands
            split_branches: if True distances and energies are split according to the
                branches. If False distances and energies are split only where branches
                are discontinuous (reducing the number of lines to plot).

        Returns:
            dict: A dictionary of the following format:
            ticks: A dict with the 'distances' at which there is a kpoint (the
            x axis) and the labels (None if no label).
            energy: A dict storing bands for spin up and spin down data
            {Spin:[np.array(nb_bands,kpoints),...]} as a list of discontinuous kpath
            of energies. The energy of multiple continuous branches are stored together.
            vbm: A list of tuples (distance,energy) marking the vbms. The
            energies are shifted with respect to the fermi level is the
            option has been selected.
            cbm: A list of tuples (distance,energy) marking the cbms. The
            energies are shifted with respect to the fermi level is the
            option has been selected.
            lattice: The reciprocal lattice.
            zero_energy: This is the energy used as zero for the plot.
            band_gap:A string indicating the band gap and its nature (empty if
            it's a metal).
            is_metal: True if the band structure is metallic (i.e., there is at
            least one band crossing the fermi level).
        """

        energies = {str(sp): [] for sp in bs.bands.keys()}

        bs_is_metal = bs.is_metal()

        if not bs_is_metal:
            vbm = bs.get_vbm()
            cbm = bs.get_cbm()

        zero_energy = 0.0
        if zero_to_efermi:
            if bs_is_metal:
                zero_energy = bs.efermi
            else:
                zero_energy = vbm["energy"]

        # rescale distances when a bs_ref is given as reference,
        # and when bs and bs_ref have different points in branches.
        # Usually bs_ref is the first one in self._bs list is bs_ref
        distances = bs.distance
        if bs_ref is not None:
            if bs_ref.branches != bs.branches:
                distances = bs._rescale_distances(bs_ref, bs)

        if split_branches:
            steps = [br["end_index"] + 1 for br in bs.branches][:-1]
        else:
            # join all the continuous branches
            # to reduce the total number of branches to plot
            steps = get_branch_steps(bs.branches)[1:-1]

        distances = np.split(distances, steps)
        for sp in bs.bands.keys():
            energies[str(sp)] = np.hsplit(bs.bands[sp] - zero_energy, steps)

        ticks = get_ticks(bs)

        vbm_plot = []
        cbm_plot = []
        bg_str = ""

        if not bs_is_metal:
            for index in cbm["kpoint_index"]:
                cbm_plot.append(
                    (
                        bs.distance[index],
                        cbm["energy"] - zero_energy if zero_to_efermi else cbm["energy"],
                    )
                )

            for index in vbm["kpoint_index"]:
                vbm_plot.append(
                    (
                        bs.distance[index],
                        vbm["energy"] - zero_energy if zero_to_efermi else vbm["energy"],
                    )
                )

            bg = bs.get_band_gap()
            direct = "Indirect"
            if bg["direct"]:
                direct = "Direct"

            bg_str = f"{direct} {bg['transition']} bandgap = {bg['energy']}"

        return {
            "ticks": ticks,
            "distances": distances,
            "energy": energies,
            "vbm": vbm_plot,
            "cbm": cbm_plot,
            "lattice": bs.lattice_rec.as_dict(),
            "zero_energy": zero_energy,
            "is_metal": bs_is_metal,
            "band_gap": bg_str,
        }

def get_branch_steps(branches):
        """
        Method to find discontinuous branches
        """
        steps = [0]
        for b1, b2 in zip(branches[:-1], branches[1:]):
            if b2["name"].split("-")[0] != b1["name"].split("-")[-1]:
                steps.append(b2["start_index"])
        steps.append(branches[-1]["end_index"] + 1)
        return steps

def get_ticks(bs):
        """
        Get all ticks and labels for a band structure plot.

        Returns:
            dict: A dictionary with 'distance': a list of distance at which
            ticks should be set and 'label': a list of label for each of those
            ticks.
        """
        bs = [bs]
        bs = bs[0] if isinstance(bs, list) else bs
        ticks, distance = [], []
        for br in bs.branches:
            s, e = br["start_index"], br["end_index"]

            labels = br["name"].split("-")

            # skip those branches with only one point
            if labels[0] == labels[1]:
                continue

            # add latex $$
            for i, l in enumerate(labels):
                if l.startswith("\\") or "_" in l:
                    labels[i] = "$" + l + "$"

            # If next branch is not continuous,
            # join the first lbl to the previous tick label
            # and add the second lbl to ticks list
            # otherwise add to ticks list both new labels.
            # Similar for distances.
            if ticks and labels[0] != ticks[-1]:
                ticks[-1] += "$\\mid$" + labels[0]
                ticks.append(labels[1])
                distance.append(bs.distance[e])
            else:
                ticks.extend(labels)
                distance.extend([bs.distance[s], bs.distance[e]])

        return {"distance": distance, "label": ticks}

def maketicks(bs, plt):
        """
        utility private method to add ticks to a band structure
        """
        ticks = get_ticks(bs)
        # Sanitize only plot the uniq values
        uniq_d = []
        uniq_l = []
        temp_ticks = list(zip(ticks["distance"], ticks["label"]))
        for i, t in enumerate(temp_ticks):
            if i == 0:
                uniq_d.append(t[0])
                uniq_l.append(t[1])
                
            else:
                uniq_d.append(t[0])
                uniq_l.append(t[1])

        plt.gca().set_xticks(uniq_d)
        plt.gca().set_xticklabels(uniq_l)

        for i in range(len(ticks["label"])):
            if ticks["label"][i] is not None:
                # don't print the same label twice
                if i != 0:
                    plt.axvline(ticks["distance"][i], color="k")
                else:
                    plt.axvline(ticks["distance"][i], color="k")
        return plt