import numpy as np
from pymatgen.core import Lattice
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import Dos, CompleteDos
from pymatgen.symmetry.kpath import KPathSetyawanCurtarolo
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

class CASTEPBands:

    def __init__(self) -> None:
        pass

    def read_file(path:str):

        def real_to_rec_lattice(real):
            reciprocal = np.zeros((3,3))
            V = np.dot(real[0],np.cross(real[1],real[2]))
            reciprocal[0] = np.cross(real[1],real[2])*(2*np.pi/V)
            reciprocal[1] = np.cross(real[2],real[0])*(2*np.pi/V)
            reciprocal[2] = np.cross(real[0],real[1])*(2*np.pi/V)
            return reciprocal;

        def get_scaled_distances(cartesian):
            moved = np.insert(cartesian, 0, cartesian[0]).reshape((np.shape(cartesian)[0]+1,np.shape(cartesian)[1]))
            cart_diff = cartesian - moved[:-1]
            distances = np.linalg.norm(cart_diff,axis=1)
            for idx, i in enumerate(distances):
                if idx == 0:
                    continue
                distances[idx] = distances[idx-1] + distances[idx]
            scaled_distances = distances/np.nanmax(distances)
            return scaled_distances;

        hartree2eV = 27.21139664
        bohr2ang = 0.529177249
        lattice, weights = [], []
        data = {}
        with open(path,'r') as f:
            lines = f.readlines()
        num_kpoints = int(lines[0].strip().split()[3])
        num_spins = int(lines[1].strip().split()[4])
        num_bands = int(lines[3].strip().split()[3])
        fermi_e = float(lines[4].strip().split()[5])*hartree2eV
        for i in range(3):
            lattice.append([float(x)*bohr2ang for x in lines[6+i].strip().split()])
        rec_lattice = real_to_rec_lattice(lattice)
        data['reciprocal_lattice'] = rec_lattice
        kpt_cart = np.zeros((num_kpoints,3))
        kpt_frac = np.zeros((num_kpoints,3))
        eigenvalues_efermi = np.zeros((num_bands, num_spins, num_kpoints))
        length_kpt_block = (num_bands+1)*num_spins+1
        for k in range(num_kpoints):
            kpt_line_index = 9+k*length_kpt_block
            frac=np.array([float(x) for x in lines[kpt_line_index].strip().split()[2:5]])
            weights.append(float(lines[kpt_line_index].strip().split()[5]))
            kpt_cart[k] = np.dot(rec_lattice,frac)
            kpt_frac[k] = frac
            for j in range(num_spins):
                for i in range(num_bands):
                    line_index = kpt_line_index + j*(num_bands+1)+i+2
                    eigenvalues_efermi[i,j,k] = float(lines[line_index].strip().split()[0])*hartree2eV-fermi_e
        data['scaled_kpt_path'] = get_scaled_distances(kpt_cart)
        data['eigenval_efermi_0'] = eigenvalues_efermi
        data['kpt_cart'] = kpt_cart
        data['kpt_frac'] = kpt_frac
        data['num_eigen'] = num_bands
        data['num_kpt'] = num_kpoints
        data['weights'] = weights
        return data;

    def read_bands2pmg(path:str=None, suffix:str=None, cell_file:str=None, export = False):

        def create_label_dict(kpath, kpts):
            naming_dict = kpath.kpath['kpoints']
            labels = {}
            high_symm_indices = []
            for key in naming_dict.keys():
                labels[key] = []
            for i in range(len(kpts)):
                for key, val in naming_dict.items():
                    if np.allclose(val,kpts[i]):
                        labels[key].append(i)
                        high_symm_indices.append(i)
            return naming_dict, high_symm_indices;

        def read_cell2pmg(path:str):
            cell = np.zeros((3,3))
            species = []
            coords = []
            with open(path,'r') as f:
                lines = f.readlines()
            
            for i,line in enumerate(lines):
                if '%BLOCK LATTICE_CART' in line or '%BLOCK lattice_cart' in line:
                    for i in range(3):
                        temp = next(f).strip().split()
                        for j in range(len(temp)):
                            cell[i,j] = float(temp[j])
                    lattice_obj = Lattice(cell)
                if '%BLOCK LATTICE_ABC' in line:
                    axes = [float(i) for i in lines[i+1].strip().split()]
                    angles = [float(i) for i in lines[i+2].strip().split()]
                    #print(axes,angles)
                    lattice_obj = Lattice.from_parameters(a= axes[0], b=axes[1], c = axes[2], alpha = angles[0], beta = angles[1], gamma = angles[2])
                if '%BLOCK POSITIONS_ABS' in line:
                    while True:
                        temp = next(f).strip().split()
                        if temp[0] == '%ENDBLOCK':
                            break
                        species.append(temp[0])
                        coords.append([float(temp[1]),float(temp[2]),float(temp[3])])
                    cartesian = True
                    break
                if '%BLOCK POSITIONS_FRAC' in line or '%BLOCK positions_frac' in line:
                    while True:
                        temp = next(f).strip().split()
                        if temp[0] == '%ENDBLOCK':
                            break
                        species.append(temp[0])
                        temp_coord = []
                        for item in temp[1:]:
                            if '/' in item:
                                fraction = item.split('/')
                                temp_coord.append(float(fraction[0])/float(fraction[1]))
                            else: temp_coord.append(float(item))
                        coords.append(temp_coord)
                    cartesian = False
                    break
            return Structure(lattice_obj,species, coords, coords_are_cartesian= cartesian);


        num_kpoints, num_spin_comps, num_electrons_up, num_electrons_down, num_bands, fermi_energy = 0,0,0,0,0,0
        kpts_coordinates = []
        eigenvalues = {}
        cell = np.zeros((3,3))
        found = False
        if path == None:
            path = f'./structures/'
        if path[-1] != '/': path += '/'
        listOfFiles = os.listdir(path)
        # create output classes for each of the output files
        for item in listOfFiles:
            # print(item.split('.')[-1])
            if not suffix == None:
                if suffix in item.split('.')[-1]:
                    bands_item = item
                    found = True
                    seed = bands_item.replace('.bands','')
            else:
                if '.bands' in item[-6:] and '.orig' not in item and '.o2b' not in item:
                    bands_item = item
                    found = True
                    seed = bands_item.replace('.bands','')
        if not found: raise EOFError('The supplied path did not contain a fitting bands file. Please ensure a proper *.bands file exists!')

        with open(f'{path}{bands_item}','r') as f:
            lines = f.readlines()

        print(f'{path}{bands_item}')
        
        for idx,line in enumerate(lines):
            line = line.strip()
            if 'Number of k-points' in line:
                num_kpoints = int(lines[idx].split()[3])
                num_spin_comps = int(lines[idx+1].split()[4])
                num_electrons = lines[idx+2].split()
                num_electrons_up = float(num_electrons[3])
                num_bands = int(lines[idx+3].split()[3])
                fermi_energy = float(lines[idx+4].split()[5])*27.2113966
                #print(fermi_energy)

                kpts_coordinates = np.zeros((num_kpoints,3))
                eigenvalues[Spin.up] = np.zeros([num_bands, num_kpoints])
                if num_spin_comps > 1:
                    num_electrons_down = float(num_electrons[4])
                    eigenvalues[Spin.down] = np.zeros([num_bands, num_kpoints])
                for i in range(3):
                    cell[i,:] = [float(x) for x in lines[idx+6+i].split()]
                print(cell)
                lattice_obj = Lattice(cell)
            if 'K-point ' in line and 'Number of k-points' not in line:
                # print('K-point') 
                temp = line.split()
                index = int(temp[1])-1
                kpts_coordinates[index] = [float(temp[2]),float(temp[3]),float(temp[4])]
                for i in range(num_bands):
                    eigenvalues[Spin.up][i][index] = float(lines[idx+2+i].strip())*27.2113966
                if num_spin_comps > 1:
                    for i in range(num_bands-1):
                        eigenvalues[Spin.down][i][index] = float([idx+2+num_bands+1+i].strip())*27.2113966

        if cell_file == None:
            for item in listOfFiles:
                if '_geometry.cell' in item:
                    kpt_path = KPathSetyawanCurtarolo(SpacegroupAnalyzer(read_cell2pmg(f'{path}{item}')).get_primitive_standard_structure())#Should use the Setyawan-Curtarolo Convention
        else: 
            kpt_path = KPathSetyawanCurtarolo(SpacegroupAnalyzer(read_cell2pmg(cell_file)).get_primitive_standard_structure())#Should use the Setyawan-Curtarolo Convention
            
        high_symm_dict, high_symm_indices = create_label_dict(kpt_path, kpts_coordinates)
        print(high_symm_dict)
        print(high_symm_indices)
        print(num_kpoints)
        final_kpt_coordinates = np.zeros((num_kpoints+len(high_symm_indices)-2,3))
        final_eigenvalues = {Spin.up : np.zeros([num_bands, num_kpoints+len(high_symm_indices)-2])}
        if num_spin_comps > 1:
            final_eigenvalues = {Spin.down : np.zeros([num_bands, num_kpoints+len(high_symm_indices)-2])}
        print(high_symm_indices)
        for i in range(len(high_symm_indices)-1):
            final_kpt_coordinates[high_symm_indices[i]+i:high_symm_indices[i+1]+1+i] = kpts_coordinates[high_symm_indices[i]:high_symm_indices[i+1]+1]
            final_eigenvalues[Spin.up][:,high_symm_indices[i]+i:high_symm_indices[i+1]+1+i] = eigenvalues[Spin.up][:,high_symm_indices[i]:high_symm_indices[i+1]+1]
            if num_spin_comps > 1:
            final_eigenvalues[Spin.down][:,high_symm_indices[i]+i:high_symm_indices[i+1]+1+i] = eigenvalues[Spin.down][:,high_symm_indices[i]:high_symm_indices[i+1]+1]
        new_bandstruct = BandStructureSymmLine(final_kpt_coordinates, final_eigenvalues,lattice_obj.reciprocal_lattice, fermi_energy,high_symm_dict,coords_are_cartesian=False);
        if export:
            with open(f'./structures/band_jsons/{seed}.json', 'w') as f:
                json.dump(new_bandstruct.as_dict(), f)
        return new_bandstruct;