from configparser import MAX_INTERPOLATION_DEPTH
from multiprocessing.spawn import old_main_modules
import warnings
class OptaDOSOutput:
    
    def __init__(self, path) -> None:
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
            if 'Electronic Data' in line:
                self.num_bands = int(lines[idx+1].strip().split()[5])
                temp = lines[idx+2].strip().split()
                if '***' not in temp:    
                    self.k_grid = [int(temp[4]),int(temp[6]),int(temp[8])]
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
            if 'Work Function' in line and not 'Effective' in line:
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
            if 'Lattice Vectors' in line:
                for i in range(3):
                    self.lattice[i] = [float(x) for x in lines[idx+1+i].strip().split()[1:]]

            if self.iprint > 1 and 'jdos_max_energy' in line:
                self.jdos_max_energy = float(line.strip().split()[3])
    
    def create_bandstructure(self,):
        
        return;