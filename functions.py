#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 10:50:04 2022

@author: brunocamino
"""
import ase

def calc_surface_energy(bulk_energy,slab_energy,n_units):
    #bulk_energy (float): energy of the bulk
    #slab_energy (float): energy of the slab
    #n_units (int): number of bulk units in the slab
    
    e_surf = 0.5*(slab_energy-(bulk_energy*n_units))
    
    return e_surf
    

def get_qe(optados_output):
    #This is for you to develop :)
    pass

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
        'orbitals2bands': '/rds/general/user/fcm19/home/modules_codes/CASTEP-19.11/obj/linux_x86_64_ifort17--serial/orbitals2bands',
        'optados': '/rds/general/user/fcm19/home/modules_codes/optados/optados.x'
    }
    
    output = ['#!/bin/bash  --login\n']
    print(output)
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
    
    with open('./{0}/{0}.qsub'.format(options['seed_name']), 'w') as f:
        for line in output:
            f.write(line)
    return;

    #based on the CASTEP calculator functionality. For documentation see 
    #https://wiki.fysik.dtu.dk/ase/ase/calculators/castep.html#module-ase.calculators.castep
def generate_castep_input(calc_struct='hello', **options): 
    #calc_struct (ase.Atoms): ASE Atoms object containing the cell and atoms for the calculation
    
    ## Optional Keywords (Must be supplied together)
    ## directory (str): Directory for generated Files
    ## label (str): seed or label for CASTEP calculation and files
    ## task (str): specification of calculation to be performed
    ## energy_cutoff (int): Energy Cutoff for Plane Waves [eV]
    ## xc_functional (str): Exchange-Correlation-Functional for the calculationa    
    ## opt_strategy (str): optimization of the calculation (default, speed, memory)
    ## spin_polarized (str): boolean value, defines if the calculation is spin_polarized, or not
    ## kpoints (str, list, dict): definition for k-point grid, has dict like functionality
    ## symmetry (str): boolean value, defines if atoms in cell should be moved to determined symmetry sites
    ## fix_all_cell (str): boolean value, defines if all the cell parameters should stay fixed during a optimisation run
    
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
        calc.param.num_dump_cycles = 0 # Prevent CASTEP from writing *wvfn* files
        
        # Define cell file options
        calc.cell.symmetry = options['symmetry']
        calc.cell.fix_all_cell = options['fix_all_cell']
        calc.set_kpts(options['kpoints'])
        
    # Prepare atoms and attach them to the current calculator
    calc_struct.calc = calc
    
    #Create input files
    calc_struct.calc.initialize()
    return;
