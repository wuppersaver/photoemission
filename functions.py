#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 10:50:04 2022

@author: brunocamino
"""
import numpy as np
import ase
import re
import warnings
from wulffpack import (SingleCrystal)
from collections import defaultdict

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
        with open('./{1}{0}/{0}.castep'.format(seed, path),'r') as f:
                 lines = f.readlines()
    else:
        with open('./structures/{0}/{0}.castep'.format(seed),'r') as f:
                 lines = f.readlines()
    for line in lines:
        if '*Warning* max. SCF cycles performed' in line:
            warnings.warn('!WARNING! The calculation in file {}.castep did not converge in the SCF cycle limit.'.format(seed))
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
    
    
    # !!Queue Configurations are specific to Cluster!!
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
    
    output = ['#!/bin/bash  --login\n']
    output.append('#PBS -N {}\n'.format(options['seed_name']))
    output.append(queues[options['queue']])
    output.append('cd $PBS_O_WORKDIR\n\nmodule load mpi intel-suite\n\n')
    output.append('PRGM={}\n\n'.format(programs[options['program']]))
    if options['bandstructure']:
        output.append('PRGMO2B={}\n\n'.format(programs['orbitals2bands']))
        output.append('CASE_O2B={}_o2b.out\n\n'.format(options['seed_name']))
    output.append('CASE_IN={0}\nCASE_OUT={0}.out\n\n'.format(options['seed_name']))
    output.append('mpiexec $PRGM $CASE_IN 2>&1 | tee $CASE_OUT\n\n')
    if options['bandstructure']:
        output.append('cp {0}.bands {0}.bands.orig\n\n'.format(options['seed_name']))
        output.append('mpiexec $PRGMO2B $CASE_IN 2>&1 | tee $CASE_O2B\n')
    
    with open('./structures/{0}/{0}.qsub'.format(options['seed_name']), 'w') as f:
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
        calc.param.cut_off_energy = str(options['energy_cutoff']) + ' eV'
        calc.param.xc_functional = options['xc_functional']
        calc.param.opt_strategy = options['opt_strategy']
        calc.param.fix_occupancy = options['fix_occup']
        calc.param.spin_polarized = options['spin_polarized']
        calc.param.max_scf_cycles = options['max_scf_cycles']
        calc.param.num_dump_cycles = 0 # Prevent CASTEP from writing *wvfn* files
        if options['continuation']:
            calc.param.continuation = 'Default'  
        if options['write_potential']:
            calc.param.write_formatted_potential = 'TRUE'
        if options['write_density']:
            calc.param.write_formatted_density = 'TRUE'
        if options['extra_bands']:
            calc.param.nextra_bands = options['number_extra_bands']
        
        # Define cell file options
        if options['snap_to_symmetry']:
            calc.cell.snap_to_symmetry = 'TRUE'
        if options['task'] == 'BandStructure':
            #print('bandstructure')
            band_path = calc_struct.cell.bandpath(options['bandstruct_path'], density = options['bandstruct_kpt_dist'])
            print(band_path)
            #print(band_path.kpts())
            calc.set_bandpath(bandpath=band_path)
            calc.cell.bs_kpoint_path_spacing = options['bandstruct_kpt_dist']
            with open(f"./structures/bandpaths/{options['label']}.bandpath", 'w') as f:
                for i in range(len(options['bandstruct_path'])):
                    if i == len(options['bandstruct_path'])-1:
                        f.write(options['bandstruct_path'][i] + '\n')
                    else:    
                        f.write(options['bandstruct_path'][i] + ',')
                f.write(np.array2string(band_path[0]))
        calc.set_kpts(options['kpoints'])
        calc.cell.fix_all_cell = options['fix_all_cell']
       
        
    # Prepare atoms and attach them to the current calculator
    calc_struct.calc = calc
    
    #Create input files
    calc_struct.calc.initialize()
    return;   
def write_bs_path(seed):
    lines = []
    path = ''
    with open('./structures/{0}/{0}.cell'.format(seed), 'r') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if 'BS_KPOINT_LIST:' in lines[i]:
            #print(lines[i])
            temp = lines[i].replace('BS_KPOINT_LIST: ', '').replace('[','').replace(']', '').replace("'", '').replace('\n','')
            path_list = list(temp.split(', '))
            print(path_list)
            path = '%BLOCK BS_KPOINT_PATH\n'
            for point in path_list:
                path += point + '\n'
            path += '%ENDBLOCK BS_KPOINT_PATH\n'
            lines[i] = path
            with open(f'./structures/{seed}/{seed}.bandpath', 'w') as f:
                print('writing_bandpath_file')
                f.write(path)
            break
    with open('./structures/{0}/{0}.cell'.format(seed), 'w') as f:
        print('writing_cell_file')
        for item in lines:
            f.write(item)
    return;
def get_wulff_fractions(mat_structure, facets_energies : dict):
    oh = SingleCrystal(facets_energies, mat_structure)
    fractions = oh.facet_fractions
    new = defaultdict(list)

    for d in (facets_energies, fractions):
        for key, value in d.items():
            new[key].append(value)
    return new;
def average_potential_from_file(input_file):

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
def create_bands_file(calc_struct,seed : str):
    calc = ase.calculators.castep.Castep(check_castep_version = False,keyword_tolerance=3)
    calc_struct.calc = calc
    bandstructure = calc_struct.calc.band_structure(bandfile=f'./structures/{seed}/{seed}.bands')
    print(bandstructure)
    bandstructure.write(fd = f'./structures/band_jsons/{seed}.json')
    return bandstructure;
def read_bands2pmg(seed : str):
    #WORK IN PROGRESS
    #Format of the pmg object:
    #classBandStructure(kpoints, eigenvals, lattice, efermi, labels_dict=None, coords_are_cartesian=False, structure=None, projections=None)
    
    import re
    import numpy as np
    from pymatgen.electronic_structure.core import Spin
    from pymatgen.core.lattice import Lattice
    from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
    import json

    lines = []
    bandpath_lines = []
    k_points_coordinates = []
    eigenvalues = {}
    labels_dict = {}
    cell = []

    # Read in files
    with open(f'./structures/{seed}/{seed}.bands','r') as f:
        lines = f.readlines()
    with open(f'./structures/bandpaths/{seed}.bandpath', 'r') as f:
        bandpath_lines = f.readlines()
    
    #Read header information
    num_kpoints = int(re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines.pop(0))[0])
    num_spin_comps = int(re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines.pop(0))[0])
    num_electrons = float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines.pop(0))[0])
    num_bands = int(re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines.pop(0))[0])
    fermi_energy = float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines.pop(0))[0]) * 27.2113966

    #Initialise the Bands and Labels dictionary
    eigenvalues[Spin.up] = np.zeros([num_bands, num_kpoints])

    #labels[Spin.up] = np.zeros([num_bands, num_kpoints])
    if num_spin_comps >1:
        eigenvalues[Spin.down] = np.zeros([num_bands, num_kpoints])

    # Read Cell vectors
    del lines[0]
    cell.append(re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines.pop(0)))
    cell.append(re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines.pop(0)))
    cell.append(re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines.pop(0)))

    # Create Pymatgen Lattice Object with Cell dimensions
    lattice_obj = Lattice(cell)

    # Read .bandpath file to create high symmetry point label dictionary
    symm_points = bandpath_lines[0].replace('G','\\Gamma').replace('\n','').split(',')
    for i in range(len(symm_points)):
        coordinates = re.findall(r"[-+]?(?:\d*\.\d+|\d+)",bandpath_lines[i+1])
        labels_dict[symm_points[i]] = [float(coordinates[0]),float(coordinates[1]),float(coordinates[0])]
    print(labels_dict)

    # Read in K-point path coordinates and Eigenvalues for all bands at each k-point
    if num_spin_comps > 1: # The ordering of K-points can be random, if more than one spin component is used.
        k_points_coordinates = np.zeros((num_kpoints,3))
        for j in range(num_kpoints):
            temp = re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines[j*(num_bands*num_spin_comps+num_spin_comps+1)])
            k_points_coordinates[int(temp[0])-1][0]=float(temp[1])
            k_points_coordinates[int(temp[0])-1][1]=float(temp[2])
            k_points_coordinates[int(temp[0])-1][2]=float(temp[3])
            for i in range(num_bands):
                eigenv_spin1_temp = re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines[j*(num_bands+num_spin_comps+1)+i+2])[0]
                eigenv_spin2_temp = re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines[j*(num_bands+num_spin_comps+1)+i+2+num_bands])[0]
                eigenvalues[Spin.up][i][int(temp[0])-1]=float(eigenv_spin1_temp)
                eigenvalues[Spin.down][i][int(temp[0])-1]=float(eigenv_spin2_temp)
    else: 
        for j in range(num_kpoints):
            temp = re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines[j*(num_bands+num_spin_comps+1)])
            k_points_coordinates.append([float(temp[1]), float(temp[2]), float(temp[3])])
            for i in range(num_bands):
                eigenv_temp = re.findall(r"[-+]?(?:\d*\.\d+|\d+)",lines[j*(num_bands+num_spin_comps+1)+i+2])[0]
                eigenvalues[Spin.up][i][j]=float(eigenv_temp)
                  
    # Create Pymatgen Bandstructure Object and save as .json file. 
    bandstruct_object = BandStructureSymmLine(kpoints = k_points_coordinates, eigenvals = eigenvalues,lattice = lattice_obj.reciprocal_lattice, efermi = fermi_energy,labels_dict = labels_dict,coords_are_cartesian=False)
    #with open(f'./structures/band_jsons/{seed}.json', 'w') as f:
    #    json.dump(bandstruct_object.as_dict(), f)       
    return bandstruct_object;