import warnings
class CastepOutput:
    def __init__(self, path) -> None:
        scf_lines = []
        self.is_converged = True

        with open(path,'r') as f:
            for line in f:
                line = line.strip()
                if '*Warning* max. SCF cycles performed' in line:
                    warnings.warn(f'!WARNING! The calculation in file {path} did not converge in the SCF cycle limit.')
                    self.is_converged = False
                if 'Final energy' in line:
                    self.ks_total_energy = float(line.split()[4])
                if 'Final free energy' in line:
                    self.mermin_free_energy = float(line.split()[5])
                if 'NB est. 0K energy' in line:
                    self.estimate_energy_0K = float(line.split()[6])
                if '-- SCF' in line:
                    scf_lines.append(line)

        self.fermi_energy = float(scf_lines[-2].split()[2])
        self.number_iterations = float(scf_lines[-2].split()[0])