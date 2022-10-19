
from ase.atoms import Atoms
from ase.build import fcc111
from pymatgen.core.structure import Structure
import warnings
class PhotoemissionCalculation:
    def __init__(self, structure:Atoms = None, general_options = None,castep_options = None,optados_options = None) -> None:
        GeneralOptions = {
            'seed' : 'Test',
            'directory': f"./structures/Test",
            'tasks' : ['Geometry','Spectral','OD_Fermi','Workfunction','OD_Photo_Sweep','BandStructure'],
            'template_scripts' : True
        }
        CastepOptions = {
            'xc_functional': 'PBE',
            'energy_cutoff': 300,
            'elec_energy_tol': 1e-8,
            'opt_strategy': 'Speed',
            'fix_occup' : False,
            'mixing_scheme' : 'Broydon', #Choose from Broydon or Pulay
            'smearing_width' : 300,
            'spectral_task' : 'All', #Choose if task is Spectral: DOS, BandStructure, Optics, CoreLoss, All
            'spin_polarized': False,
            'max_scf_cycles': 1000,
            'write_potential': False,
            'write_density': False,
            'extra_bands': False,
            #Cell File Instructions
            'kpoints': (9,9,9),
            'generate_symmetry' : False,
            'snap_to_symmetry': False,
            'fix_all_cell': False,
            'continuation': False,
            'bandstruct_path': 'GM',
            'bandstruct_kpt_dist': 0.0184,
            'spectral_kpt_grid': (9,9,9),
            'calculate_pdos': False 
        }
        OptaDosOptions = {
            'optados_task': 'pdos', # Choose: dos(default), compare_dos, compare_jdos, jdos, pdos, optics, core, all
            'broadening': 'adaptive', #Choose: adaptive(default), fixed, linear
            'iprint': '1', #Choose: 1 - bare minimum, 2 - with progress reports, 3 - fulld debug output
            'efermi': 'optados', #Choose: optados - recalculate, file - read from CASTEP file, insulator - count filled bands, float - supplied by user
            'dos_spacing': '0.001', #DOS spacing in unit (default: eV): default - 0.1
            'pdos': 'angular', #Choose: angular, species_ang, species, sites or more detailed descriptions such as: 
            #PDOS : sum:Si1-2(s) - sum of s-chnnls on 2 Si atms (1 proj), 
            #PDOS : Si1;Si2(s) - DOS on Si atom 1 and DOS on s-channel of Si atom 2 (2 proj) 
            'photo_options': {
                'work_function' : 4.556,
                'surface_area' : 6.339744873,
                'slab_volume' : 190.942338,
                'elec_field' : 0,
                'imfp_const' : 19.0,
                'JDOS_SPACING' : 0.1,
                'JDOS_MAX_ENERGY' : 25,
                'BROADENING' : 'linear',
                'OPTICS_GEOM' : 'unpolar',
                'optics_qdir' : [1, 1.000, 1.0000],
                'photon_energy' : 21.2,
                'linear_smearing' : 0.026,
                'fixed_smearing' :  0.026,
                'optics_intraband' : True,
                'photo_model' : '1step',
                'momentum' : 'crystal',
                'hybrid_linear' : True,
                'temp' : 300,
                'theta_lower' : 59,
                'theta_upper' : 61,
                'iprint' : 1,
                },
            'sweep_options': {
                'parameter' : 'photon_energy', # choose photon_energy or temp or elec_field or work_function
                'min' : 4,
                'stepsize' : 0.1,
                'max' : 5.5,
                },
        }
        if structure == None : structure = fcc111('Cu', size=(1,1,5), vacuum=20.0)
        if general_options == None : general_options = GeneralOptions
        if castep_options == None : castep_options = CastepOptions
        if optados_options == None : optados_options = OptaDosOptions
        self.options = {'general':general_options,'structure':structure,'castep':castep_options,'optados':optados_options}
