import re
# from configparser import MAX_INTERPOLATION_DEPTH
# from multiprocessing.spawn import old_main_modules
import warnings
class OptaDOSOutput:
    
    def __init__(self, path:str) -> None:
        scf_lines = []
        self.path = path
        #self.seed = path.split('/')[-1].split('.')[0]
        self.lattice = [[],[],[]]
        self.atoms,self.atom_coordinates = [],[]
        self.layers= {}
        optados_jobs = ['output_DOS', 'output_partial_DOS', 'output_projected_BS', 'output_jDOS', 'output_optic_response', 'output_corelvl_spec', 'photoemission']
        self.methods = []
        lines = None
        with open(path,'r') as f:
            lines = f.readlines()            
        for idx,line in enumerate(lines):
            line = line.strip()
            if 'JOB CONTROL' in line:
                for i in range(7):
                    if 'True' in lines[idx+1+i]:
                        self.methods.append(optados_jobs[i])
                #print(next(f).strip().split()[4])
                self.iprint = int(lines[idx+8].strip().split()[4])
                if 'True' in lines[idx+9]: self.file_format = 'old'
                else: self.file_format = 'new'
            if 'Smearing Width' in line:
                self.temperature = float(line.split()[4])*11604.51812
            
            if 'PHOTOEMISSION PARAMETERS' in line:
                temp_lines = lines[idx+1:idx+16]
                for line in temp_lines:
                    if 'Photoemission Model' in line: self.photo_model = ' '.join(line.split()[4:5])
                    if 'Photoemission Final State' in line: self.photo_final_state = ' '.join(line.split()[5:7])
                    if 'Photon Energy' in line : self.photo_photon_energy = float(line.split()[5])
                    if 'Work Function' in line : self.photo_work_function = float(line.split()[5])
                    if 'Surface Area'  in line : self.photo_surface_area = float(line.split()[5])
                    if 'Slab Volume'   in line : self.photo_slab_volume = float(line.split()[5])
                    if 'IMFP Constant' in line : self.photo_imfp = float(line.split()[5])
                    if 'Electric Field Strength' in line : self.photo_elec_field = float(line.split()[6])
                    if 'Smearing Temperature' in line : self.photo_smearing_temperature = float(line.split()[5])
            
            if 'Lattice Vectors' in line:
                for i in range(3):
                    self.lattice[i] = [float(x) for x in lines[idx+1+i].strip().split()[1:]]

            if 'Electronic Data' in line:
                self.num_bands = int(lines[idx+1].strip().split()[5])
                temp = lines[idx+2].strip().split()
                if '***' not in temp:
                    for index,item in enumerate(temp):
                        temp[index] = item.replace('x','')
                    if int(temp[4]) > 99:  self.k_grid = [int(temp[4]),int(temp[5]),int(temp[7])]
                    else: self.k_grid = [int(temp[4]),int(temp[6]),int(temp[8])]
                self.num_kpts = int(lines[idx+3].strip().split()[5])
                if 'True' in lines[idx+4].strip(): self.spin_polarised = True
                else: self.spin_polarised = False
                self.num_electrons = float(lines[idx+5].strip().split()[5])

            if 'Projected Density Of States Calculation' in line:
                #print('fermi analysis')
                for i in range(6):
                    line = next(f)
                line = lines[idx+6].strip().split()
                #print(line)
                self.fermi_e = float(line[6])

            if 'Max number of atoms' in line:
                self.max_atoms = int(line.split()[5])
                self.max_layers = int(line.split()[10])
            if 'Work Function' in line and 'Photon Energy' in line:
                self.work_fct = float(line.split()[3])
                self.photon_e = float(line.split()[7])
                temp = lines[idx+1].strip().split()
                #print(temp)
                self.work_fct_effect = float(temp[4])
                self.e_field = float(temp[8].replace('V/A',''))
                temp = lines[idx+1]
                if 'Free electron state' in temp: self.final_state = 'free electron'
                if 'Bloch state' in temp: self.final_state = 'bloch'
                #print(self.final_state)
            if '  Quantum Efficiency  ' in line:
                #print(self.max_layers)
                for i in range(self.max_layers):
                    temp = lines[idx+1+i].strip().split()
                    self.layers[str(i)] = [temp[1], int(temp[2]), int(temp[3])]
            if 'Total Quantum Efficiency' in line:
                self.qe = float(line.split()[5])
                self.mte = float(lines[idx+1].split()[6])
            #print(line)
            if hasattr(self, 'iprint'):
                if self.iprint > 1 and 'jdos_max_energy' in line:
                    self.jdos_max_energy = float(line.strip().split()[3])
        #print ([value for value in self.__dict__.keys()])

    def create_bandstructure(self,):
        
        return;