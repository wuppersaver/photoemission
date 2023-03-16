import os
import numpy as np
import pandas as pd
from pymatgen.electronic_structure.dos import *

class PDosDataFrame:
    def __init__(self, path:str, choice:str, spec:str=None) -> None:
        choices = {
            'site' : 'sites',
            'orbital': 'angular',
        }
        if choice == 'specified' and spec != None:
            choices['specified'] = spec
        energies, total, species = [],[],[]
        columns, channels, column_keys, totals = [], [], {}, {Spin.up:[], Spin.down:[]}
        header_string = '#+----------------------------------------------------------------------------+'
        header, values = [],[]
        spin_channels = False
        shifted_efermi = False
        self.path = path
        listOfFiles = os.listdir(path)
         # create output classes for each of the output files
        for item in listOfFiles:
            if '_fermi.odo' in item:
                with open(path + item,'r') as g:
                    for line in g:
                        if 'Shift energy scale so fermi_energy=0' in line: shifted_efermi = bool(line.split()[7]=='True')
                        if 'Projected Density Of States Calculation' in line:
                            for i in range(3): line = next(g)
                            while 'Fermi energy (' not in line: 
                                line = next(g)    
                            efermi = float(line.split()[6])
                            #print('No fermi energy found, check cursor position!')
            if f'{choices[choice]}.pdos.dat' in item:
                values = np.loadtxt(path+item)
                # print(values)
                with open(path+item,'r') as f:
                    for line in f:
                        line_values = line.strip().split()
                        if '#' in line_values[0]:
                            header.append(line_values)
                        else: break
        for i, item in enumerate(header):
            if 'Column:' in item:
                column = {}
                species_temp = {}
                if 'Spin' in header[i+1]: 
                    spin_channels = True
                    print(spin_channels)
                t = i+2
                while header_string not in header[t]:
                    atom = ''.join(header[t][1:3])
                    specie = header[t][1]
                    if specie not in species_temp: species_temp[f'{specie}'] = []
                    if atom in column:
                        column[atom].append(header[t][3])
                        if spin_channels: channels.append(header[t][4])
                    else: 
                        column[atom] = []
                        column[atom].append(header[t][3])
                        if spin_channels: channels.append(header[t][4])
                    t+=1
                temp = ''
                for atom in column.keys():
                    temp += f'{atom},'
                    for item in column[atom]:
                        temp += item
                    temp += ' '
                columns.append(temp[:-1])
                temp = ''
                for specie in species_temp.keys():
                    temp += f'{specie}'
                species.append(temp)
        data = []
        self.total_dos = []
        sums = sum(values[:,1:])
        normed = np.zeros(values.shape)
        for idx,summ in enumerate(sums):
            if summ > 1E-100:
                normed[:,idx] = values[:,idx]/np.nanmax(values[:,idx],axis=0)
        normed[:,0] = values[:,0]
        data_norm = []
        for row in values:
            if sum(abs(row[1:]))>0:
                if shifted_efermi: 
                    self.total_dos.append([round(row[0],3),sum(abs(row[1:]))])
                else: self.total_dos.append([round(row[0]-efermi,3),sum(abs(row[1:]))])
                for idx,item in enumerate(row[1:]):
                    if shifted_efermi: 
                        if spin_channels: 
                            data.append([species[idx],columns[idx],channels[idx],round(row[0],3),item])
                        else: data.append([species[idx],columns[idx],round(row[0],3),item])
                    else: 
                        if spin_channels: data.append([species[idx],columns[idx],channels[idx],round(row[0]-efermi,3),item])
                        else: data.append([species[idx],columns[idx],round(row[0]-efermi,3),item])
        for row in normed:
            if sum(abs(row[1:]))>0:
                for idx,item in enumerate(row[1:]):
                    if shifted_efermi: 
                        if spin_channels: 
                            data_norm.append([species[idx],columns[idx],channels[idx],row[0],item])
                        else: data_norm.append([species[idx],columns[idx],row[0],item])
                    else: 
                        if spin_channels: data_norm.append([species[idx],columns[idx-1],channels[idx],row[0]-efermi,item])
                        else: data_norm.append([species[idx],columns[idx],row[0]-efermi,item])
        if spin_channels: 
            self.df = pd.DataFrame(data,columns = ['Species','Site','Spin','Energy[eV]','pDOS'])  
            self.df_norm = pd.DataFrame(data_norm,columns = ['Species','Site','Spin','Energy[eV]','pDOS(normalised)'])
        else: 
            self.df = pd.DataFrame(data,columns = ['Species','Site','Energy[eV]','pDOS'])
            self.df_norm = pd.DataFrame(data_norm,columns = ['Species','Site','Energy[eV]','pDOS(normalised)'])
        self.total_dos = np.array(self.total_dos)
        self.total_dos_norm =  self.total_dos / np.nanmax(self.total_dos[:,1])
        self.total_dos_norm[:,0] = self.total_dos[:,0]
        return;