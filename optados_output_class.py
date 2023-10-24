import re
# from configparser import MAX_INTERPOLATION_DEPTH
# from multiprocessing.spawn import old_main_modules
import warnings
import numpy as np
from pathlib import Path

class OptaDOSOutput:
    
    def __init__(self, path:Path) -> None:        
        total_qe = re.compile('Total Quantum Efficiency \(electrons/photon\):\s+([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))?\s+')
        mte = re.compile('Weighted Mean Transverse Energy \(eV\):\s+([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))?\s+')
        bulk =  re.compile('Bulk\s+([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))? ')
        beginning_qe = re.compile('\| Atom \|\s+Atom Order\s+\|\s+Layer\s+\|\s+Quantum Efficiency\s+\|')
        workf_photoe = re.compile('Work Function\s+([0-9]*\.[0-9]+) eV\s+Photon Energy\s+([0-9]*\.[0-9]+) eV')
        effworkf_efield = re.compile('Effective Work Function\s+([0-9]*\.[0-9]+)\seV\s+Electric Field\s+([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))? V/A\s+')
        final_state = re.compile("Final State : ([A-Za-z]+( [A-Za-z]+)+)")
        kpts_grid = re.compile("\s+Grid size\s+:\s+([0-9]+)\s*x\s*([0-9]+)\s*x\s*([0-9]+)\s+")
        optados_jobs = ['output_DOS', 'output_partial_DOS', 'output_projected_BS', 'output_jDOS', 'output_optic_response', 'output_corelvl_spec', 'photoemission']
        self.od_parameters = {
            'iprint' : None,
            'file_format' : None,
            'smearing_width' : None,
            'optics_geom' : None,
            'q_vector' : None,
            'intraband' : None,
            'drude_broad' : None,
            'photon_energies' : [],
            'file_format' : 'new',
            'methods' :[],
        }
        self.system = {
            'lattice' : [[],[],[]],
            'atoms' : [], 
            'atom_coordinates' : [],
            'max_atoms' : 1,
            'max_layers' : 1
        }
        self.path = path
        self.qe_data = {}
        self.data_min = []
        lines = None
        with path.open() as f:
            lines = f.readlines()            
        for idx,line in enumerate(lines):
            line = line.strip()
            if 'JOB CONTROL' in line:
                for i in range(7):
                    if 'True' in lines[idx+1+i]:
                        self.od_parameters['methods'].append(optados_jobs[i])
                self.od_parameters['iprint'] = int(lines[idx+8].strip().split()[4])
                if 'True' in lines[idx+9]: self.od_parameters['file_format'] = 'old'
            if 'Smearing Width' in line:
                self.od_parameters['smearing_width'] = float(line.split()[4])
            
            if 'OPTICS' in line:
                temp_lines = lines[idx+1:idx+15]
                for line in temp_lines:
                    if 'Geometry for' in line: self.od_parameters['optics_geom'] = line.split()[6].lower()
                    if 'qvector' in line: self.od_parameters['q_vector'] = [float(x) for x in line.split()[6:9]]
                    if 'Intraband Contribution' in line: self.od_parameters['intraband'] = True
                    if 'Drude Broadening' in line: self.od_parameters['drude_broad'] = float(line.split()[4])
            if 'PHOTOEMISSION PARAMETERS' in line:
                temp_lines = lines[idx+1:idx+16]
                for line in temp_lines:
                    match = re.search('Photoemission Model\s+:\s+([A-Za-z0-9]+-Step Model)\s+',line)
                    if match:self.od_parameters['photo_model'] = match.groups()[0].lower()
                    match = re.search('Photoemission Final State\s+:\s+([A-Za-z0-9]+( [A-Za-z0-9]+)+)',line)
                    if match:self.od_parameters['photo_final_state'] = match.groups()[0]
                    match = re.search('Photon Energy\s+:\s+([0-9]*\.[0-9]+)\s+eV\s ',line)
                    if match: 
                        self.od_parameters['photoemission'] = 'single'
                        self.od_parameters['photo_photon_energy'] = float(match.groups()[0])
                    match = re.search('Photon Energy Sweep\s+:\s+([0-9]*\.[0-9]+)\s+->\s+([0-9]*\.[0-9]+) eV\s ',line)
                    if match:
                        self.od_parameters['photoemission'] = 'sweep'
                        self.od_parameters['photo_photon_min'] = float(match.groups()[0])
                        self.od_parameters['photo_photon_max'] = float(match.groups()[1])
                    match = re.search('Work Function\s+\(eV\)\s+:\s+([0-9]*\.[0-9]+)\s+',line)
                    if match : self.od_parameters['photo_work_function'] = float(match.groups()[0])
                    match = re.search('Surface Area\s+\(Ang\*\*2\)\s+:\s+([0-9]*\.[0-9]+)\s+',line)
                    if match :  self.od_parameters['photo_surface_area'] = float(match.groups()[0])
                    match = re.search('Slab Volume\s+\(Ang\*\*3\)\s+:\s+([0-9]*\.[0-9]+)\s+',line)
                    if match :  self.od_parameters['photo_slab_volume'] = float(match.groups()[0])
                    match = re.search('IMFP Constant\s+\(Ang\)\s+:\s+([0-9]*\.[0-9]+)\s+',line)
                    if match :  self.od_parameters['photo_imfp'] = float(match.groups()[0])
                    match = re.search('Electric Field Strength\s+\(V/Ang\)\s+:\s+([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))?',line)
                    if match :  self.od_parameters['photo_elec_field'] = float(match.groups()[0])
                    match = re.search('Smearing Temperature\s+\(K\)\s+:\s+([0-9]*\.[0-9]+)\s+',line)
                    if match :  self.od_parameters['photo_smearing_temperature'] = float(match.groups()[0])
            
            if 'Lattice Vectors' in line:
                for i in range(3):
                    self.system['lattice'][i] = [float(x) for x in lines[idx+1+i].strip().split()[1:]]

            if 'Electronic Data' in line:
                self.system['num_bands'] = int(lines[idx+1].strip().split()[5])
                temp = lines[idx+2].strip()
                if '***' not in temp:
                    self.system['k_grid'] = [int(x) for x in kpts_grid.search(temp).groups()]
                self.system['num_kpts'] = int(lines[idx+3].strip().split()[5])
                if 'True' in lines[idx+4].strip(): self.system['spin_polarised'] = True
                else: self.system['spin_polarised'] = False
                self.system['num_electrons'] = float(lines[idx+5].strip().split()[5])

            if 'Projected Density Of States Calculation' in line:
                #print('fermi analysis')
                # for i in range(6):
                #     line = next(f)
                line = lines[idx+6].strip().split()
                self.system['fermi_e'] = float(line[6])

            if 'Max number of atoms' in line:
                self.system['max_atoms'] = int(line.split()[5])
                self.system['max_layers'] = int(line.split()[10])

            match = workf_photoe.search(line)
            if match:
                self.od_parameters['workfct'] = float(match.groups()[0])
                photon_energy = match.groups()[1]
                self.qe_data[photon_energy] = {'layers'  : np.zeros((self.system['max_layers'])),
                                          'atoms'   : [],
                                          'bulk'    : 0,
                                          'total'   : 0,
                                          'mte'     : 0 }
                self.od_parameters['photon_energies'].append(float(match.groups()[1]))
            match = effworkf_efield.search(line)
            if match:
                # print(match.groups())
                self.od_parameters['elec_field_strength'] = match.groups()
            match = beginning_qe.match(line)
            if match:
                qe_temp = lines[idx+1:idx+1+self.system['max_atoms']]
                for atom in qe_temp:
                    info = atom.split()[1:-1]
                    info[1] = int(info[1])
                    info[2] = int(info[2])
                    info[3] = float(info[3])
                    self.qe_data[photon_energy]['layers'][info[2]-1] += info[3]
                    self.qe_data[photon_energy]['atoms'].append(info)
            match = bulk.search(line)
            if match:
                self.qe_data[photon_energy]['bulk'] = float('E'.join(match.groups()))
            match = total_qe.search(line)
            if match:
                self.qe_data[photon_energy]['total'] = float('E'.join(match.groups()))
                self.data_min.append([float(photon_energy),float('E'.join(match.groups())),0])
            match = mte.search(line)
            if match:
                self.qe_data[photon_energy]['mte'] = float('E'.join(match.groups()))
                self.data_min[-1][-1] = float('E'.join(match.groups()))
            if hasattr(self, 'iprint'):
                if self.iprint > 1 and 'jdos_max_energy' in line:
                    self.jdos_max_energy = float(line.strip().split()[3])
        self.data_min = np.array(self.data_min)
    def create_bandstructure(self,):
        
        return;