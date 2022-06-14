from configparser import MAX_INTERPOLATION_DEPTH
from multiprocessing.spawn import old_main_modules
import warnings
class OptaDOSOutput:
    
    def __init__(self, path) -> None:
        scf_lines = []
        self.path = path
        #self.seed = path.split('/')[-1].split('.')[0]
        lattice = [[],[],[]]
        atoms,atom_coordinates = [],[]
        self.layers= {}
        optados_jobs = ['output_DOS', 'output_partial_DOS', 'output_projected_BS', 'output_jDOS', 'output_optic_response', 'output_corelvl_spec', 'photoemission']
        self.methods = []
        with open(path,'r') as f:
            for line in f:
                line = line.strip()
                if 'JOB CONTROL' in line:
                    for i in range(7):
                        line = next(f)
                        #print(line)
                        if 'True' in line:
                            self.methods.append(optados_jobs[i])
                    #print(next(f).strip().split()[4])
                    self.iprint = int(next(f).strip().split()[4])
                    if 'True' in next(f): self.file_format = 'old'
                    else: self.file_format = 'new'
                if 'Smearing Width' in line:
                    self.temperature = float(line.split()[4])*11604.51812
                if 'Electronic Data' in line:
                    self.num_bands = int(next(f).strip().split()[5])
                    temp = next(f).strip().split()
                    self.k_grid = [int(temp[4]),int(temp[6]),int(temp[8])]
                    self.num_kpts = int(next(f).strip().split()[5])
                    if 'True' in next(f).strip(): self.spin_polarised = True
                    else: self.spin_polarised = False
                    self.num_electrons = float(next(f).strip().split()[5])

                if 'Projected Density Of States Calculation' in line:
                    #print('fermi analysis')
                    for i in range(6):
                        line = next(f)
                    line = line.strip().split()
                    #print(line)
                    self.fermi_e = float(line[6])

                if 'Max number of atoms' in line:
                    self.max_atoms = int(line.split()[5])
                    self.max_layers = int(line.split()[10])
                if 'Work Function' in line:
                    self.work_fct = float(line.split()[3])
                    self.photon_e = float(line.split()[7])
                    temp = next(f).strip().split()
                    self.work_fct_effect = float(temp[4])
                    self.e_field = float(temp[8].replace('V/A',''))
                    temp = next(f)
                    if 'Free electron state' in temp: self.final_state = 'free electron'
                    if 'Bloch state' in temp: self.final_state = 'bloch'
                    #print(self.final_state)
                if 'Quantum Efficiency' in line:
                    #print(self.max_layers)
                    for i in range(self.max_layers):
                        temp = next(f).strip().split()
                        self.layers[str(i)] = [temp[1], int(temp[2]), int(temp[3])]
                if 'Total quantum efficiency' in line:
                    self.qe = float(line.split()[5])
                    temp = next(f).strip()
                    if 'NaN' in temp: self.mte = 0
                    else: self.mte = float(temp.split()[6])
                if 'Lattice Vectors' in line:
                    for i in range(3):
                        temp = next(f).strip().split()
                        lattice[i] = [float(x) for x in temp[1:4]]
    
    def create_bandstructure(self,):
        
        return;