from pymatgen.core.structure import Structure
import warnings
class CastepOutput:
    
    def __init__(self, path) -> None:
        scf_lines = []
        self.path = path
        self.seed = path.split('/')[-1].split('.')[0]
        lattice = [[],[],[]]
        atoms,atom_coordinates = [],[]
        with open(path,'r') as f:
            for line in f:
                line = line.strip()
                if '*Warning* max. SCF cycles performed' in line:
                    warnings.warn(f'!WARNING! The calculation in file {path} did not converge in the SCF cycle limit.')
                    self.is_converged = False
                else: self.is_converged = True
                if 'Final energy' in line:
                    self.ks_total_energy = float(line.split()[4])
                if 'Final free energy' in line:
                    self.mermin_free_energy = float(line.split()[5])
                if 'NB est. 0K energy' in line:
                    self.estimate_energy_0K = float(line.split()[6])
                if '-- SCF' in line:
                    scf_lines.append(line)
                if 'Real Lattice' in line:
                    for i in range(3):
                        temp = next(f).strip().split()
                        lattice[i] = [float(x) for x in temp[0:3]]
                if 'Total number of ions in cell' in line:
                    num_atoms = int(line.split()[7])
                if 'Fractional coordinates of atoms' in line:
                    next(f)
                    next(f)
                    for i in range(num_atoms):
                        line = next(f).strip().split()
                        atoms.append(line[1])
                        atom_coordinates.append([float(x) for x in line[3:6]])
        self.structure = Structure(lattice, atoms, atom_coordinates, coords_are_cartesian= False)       
        self.fermi_energy = float(scf_lines[-2].split()[2])
        self.number_iterations = float(scf_lines[-2].split()[0])
    
    def create_bandstructure(self,):
        
        return;