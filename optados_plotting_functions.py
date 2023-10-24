import mmap
import numpy as np

def read_in_debug_file(filename):
    mat_shape = []
    with open(filename,'r') as f:
        next(f)
        line= next(f)
        mat_shape = [int(x) for x in line.split()]
        print("Will return matrix of shape:", mat_shape)
    matrix = np.genfromtxt(filename,skip_header=2,delimiter=' ')
    matrix = np.reshape(matrix,mat_shape, order='F')
    return matrix;
    
def read_qe_calculation_values(filename:str, model:str):
    data = {}
    with open(filename, "r+b") as input:
        mm = mmap.mmap(input.fileno(),0)
        iterator = iter(mm.readline, b"")
        for line in iterator:
            line_temp = line.decode("ANSI").strip()
            if "- Printing list of values going into" in line_temp:
                data['printed_quantities'] = next(iterator).decode("ANSI").strip().split(' - ')
                if model == '3step':
                    data['e_fermi'] = float(next(iterator).decode("ANSI").strip().split()[2])
                    data['data_shape'] = [int(x) for x in next(iterator).decode("ANSI").strip().split()[2:-1]]
                    # print(data['data_shape'])
                    # data['data_shape'] = data['data_shape'][1:] + [data['data_shape'][0]]
                    # print(data['data_shape'])
                if model == '1step':
                    data['data_shape'] = [int(x) for x in next(iterator).decode("ANSI").strip().split()[2:]]
                    data['data_shape'] = data['data_shape'][1:] + [data['data_shape'][0]]
                    # print(data['data_shape'])
                data['values'],data['indices'],data['qe_values'] = [],[],[]
                line_temp = next(iterator).decode("ANSI").strip()
                while "- Finished Printing -" not in line_temp:
                    data['indices'].append([int(x) for x in line_temp.split()[1:-1]])
                    line_temp = next(iterator).decode("ANSI").strip()
                    data['qe_values'].append(float(line_temp.split()[0]))
                    data['values'].append([float(x) for x in line_temp.split()[1:]])
                    line_temp = next(iterator).decode("ANSI").strip()
        mm.close()

        data['values'] = np.array(data['values'])
        data['indices'] = np.array(data['indices'])
        data['qe_values'] = np.array(data['qe_values'])
        #print(data['values'].shape)
        #print(data['values'])
        #data['values'] = np.reshape(data['values'],data['data_shape'])
    return data;
    
# def filter_3step_values(data_in:dict):
#     '''This function takes in a dictionary of raw data and filters out the bands, where the final state is below the\n
#     fermi level and thus assumed to be occupied. This selection rule is a WIP. The information for the atom,\n
#     final state, initial state, spin, and k point are preserved and returned as dict['indices'].'''
#     data_out = {}
#     data_out['printed_quantities'] = data_in['printed_quantities']
#     data_out['original_shape'] = data_in['data_shape']
#     indices = np.empty(data_in['data_shape'][:-1],dtype=object)
#     for m in range( data_in['data_shape'][4]): # K point
#         for l in range( data_in['data_shape'][3]): # Spin
#             for k in range( data_in['data_shape'][2]): # Final State
#                 for j in range( data_in['data_shape'][1]): # Initial State
#                     for i in range( data_in['data_shape'][0]): # Layer
#                         indices[i,j,k,l,m] = f"{i}-{j}-{k}-{l}-{m}"
#     indices = indices.flatten(order='F')
#     filter_array = data_in['values'][:,0] > data_in['e_fermi']
#     data_out['values'] = data_in['values'][filter_array]
#     data_out['indices'] = indices[filter_array]
#     return data_out;

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
    
def read_bands_file(path:str):
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

def read_omes(filename:str):
    data = {}
    with open(filename, "r+b") as input:
            mm = mmap.mmap(input.fileno(),0)
            iterator = iter(mm.readline, b"")
            for line in iterator:
                line_temp = line.decode("ANSI").strip()
                if 'Final Optical Matrix Elements' in line_temp:
                    line_temp = next(iterator).decode("ANSI")
                    data['kpoint_coords'] = []
                    data['omes'] = []
                    while not 'Calculation re-parallelised over' in line_temp:
                        data['kpoint_coords'].append([float(x) for x in line_temp.strip().split()[1:]])
                        num_eigen =  int(next(iterator).decode("ANSI").strip().split()[3])
                        ome_ktemp = []
                        for index in range(num_eigen):
                            next(iterator)
                            next(iterator)
                            elements = [float(x) for x in next(iterator).decode("ANSI").strip().split()]
                            x = complex(elements[0],elements[1])
                            y = complex(elements[2],elements[3])
                            z = complex(elements[4],elements[5])
                            abs_temp = np.sqrt(float(abs(x) + abs(y) + abs(z)))
                            ome_ktemp.append([x,y,z,abs_temp])
                        data['omes'].append(ome_ktemp)
                        line_temp = next(iterator).decode("ANSI")
            mm.close()
            data['kpoint_coords'] = np.array(data['kpoint_coords'])
            data['omes'] = np.array(data['omes'])
    return data;

def read_qe_matrix_values_odo(filename:str):
    data = {}
    with open(filename, "r+b") as input:
            mm = mmap.mmap(input.fileno(),0)
            iterator = iter(mm.readline, b"")
            for line in iterator:
                line_temp = line.decode("ANSI").strip()
                if 'K-Points in Cartesian Coordinates' in line_temp:
                    line = next(iterator).decode("ANSI")
                    kpoints = []
                    while not 'Finished Printing' in line:
                        #print(line)
                        kpoints.append([float(x) for x in line.strip().split()])
                        line = next(iterator).decode("ANSI")
                    data['kpoints'] = np.array(kpoints)
                if 'QE Matrix --' in line_temp:
                    data['model'] = line_temp.strip().split()[2]
                    temp = next(iterator)
                    shape = [int(x) for x in next(iterator).decode("ANSI").strip().split()]
                    matrix = []
                    line = next(iterator).decode("ANSI")
                    while not 'Finished Printing' in line:
                        matrix.append([float(x) for x in line.strip().split()]) 
                        line = next(iterator).decode("ANSI")
                    matrix = np.array(matrix)
            mm.close()
            data['matrix'] = np.reshape(matrix.flatten(order='C'),shape,order = 'F')
    return data;


def mirror_xy(array_in):
    indices = []
    bands_added = np.zeros((1,2),dtype=float)
    for index,item in enumerate(array_in):
        reversed_x = [[-1*item[0],item[1]]]
        reversed_y = [[item[0],-1*item[1]]]
        #both = [[-1*item[0],-1*item[1]]]
        if not (array_in[:, None] == reversed_x[0]).all(-1).any() and not (bands_added[:, None] == reversed_x[0]).all(-1).any(): 
            bands_added = np.append(bands_added,reversed_x,axis=0)
            indices.append([index]+reversed_x)
        if not (array_in[:, None] == reversed_y[0]).all(-1).any() and not (bands_added[:, None] == reversed_y[0]).all(-1).any(): 
            bands_added = np.append(bands_added,reversed_y,axis=0)
            indices.append([index]+reversed_y)
        # if not (array_in[:, None] == both[0]).all(-1).any() and not (bands_added[:, None] == both[0]).all(-1).any(): 
        #     bands_added = np.append(bands_added,reversed_y,axis=0)
        #     indices.append(index)
    sort_bands = np.full(bands_added.shape,True)
    sort_bands[0,:] = [False,False]
    length=len(bands_added)
    bands_added = bands_added[sort_bands]
    bands_added = bands_added.reshape(length-1,2)
    return np.vstack((array_in,bands_added)),indices
    
def copy_mirrored_data(data_in,indices):
    shape_in = data_in.shape
    data_out = np.zeros((shape_in[0],shape_in[1]+len(indices),shape_in[2]))
    data_out[:,:shape_in[1],:] = data_in
    for idx, item in enumerate(indices):
        data_out[:,shape_in[1]+idx,:] = data_in[:,item[0],:]
        data_out[:,shape_in[1]+idx,:2] = item[1]
    return data_out;
    