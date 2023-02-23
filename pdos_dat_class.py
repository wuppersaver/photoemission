import os
import numpy as np
import pandas as pd
from pymatgen.electronic_structure.dos import *

class PDosDataFrame:
    def __init__(self, path:str) -> None:
        energies, total= [],[]
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
            if '.pdos.dat' in item:
                values = np.loadtxt(path+item)
                with open(path+item,'r') as f:
                    for line in f:
                        line_values = line.strip().split()
                        if '#' in line_values[0]:
                            header.append(line_values)
                        else: break
        for i, item in enumerate(header):
            if 'Column:' in item:
                column = {}
                if 'Spin' in header[i+1]: 
                    spin_channels = True
                    print(spin_channels)
                t = i+2
                while header_string not in header[t]:
                    atom = ''.join(header[t][1:3])
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
                columns.append(temp)
        data = []
        self.total_dos = []
        normed = values/np.nanmax(values,axis=0)
        normed[:,0] = values[:,0]
        data_norm = []
        for row in values:
            if sum(abs(row[1:]))>0:
                if shifted_efermi: 
                    self.total_dos.append([row[0],sum(abs(row[1:]))])
                else: self.total_dos.append([row[0]-efermi,sum(abs(row[1:]))])
                for idx,item in enumerate(row[1:]):
                    if shifted_efermi: 
                        if spin_channels: 
                            data.append([columns[idx],channels[idx],row[0],item])
                        else: data.append([columns[idx],row[0],item])
                    else: 
                        if spin_channels: data.append([columns[idx],channels[idx],row[0]-efermi,item])
                        else: data.append([columns[idx],row[0]-efermi,item])
        for row in normed:
            if sum(abs(row[1:]))>0:
                for idx,item in enumerate(row[1:]):
                    if shifted_efermi: 
                        if spin_channels: 
                            data_norm.append([columns[idx],channels[idx],row[0],item])
                        else: data_norm.append([columns[idx],row[0],item])
                    else: 
                        if spin_channels: data_norm.append([columns[idx-1],channels[idx],row[0]-efermi,item])
                        else: data_norm.append([columns[idx],row[0]-efermi,item])
        if spin_channels: 
            self.df = pd.DataFrame(data,columns = ['Site','Spin','Energy[eV]','pDOS'])  
            self.df_norm = pd.DataFrame(data_norm,columns = ['Site','Spin','Energy[eV]','pDOS(normalised)'])
        else: 
            self.df = pd.DataFrame(data,columns = ['Site','Energy[eV]','pDOS'])
            self.df_norm = pd.DataFrame(data_norm,columns = ['Site','Energy[eV]','pDOS(normalised)'])
        self.total_dos = np.array(self.total_dos)
        self.total_dos_norm =  self.total_dos / np.nanmax(self.total_dos[:,1])
        self.total_dos_norm[:,0] = self.total_dos[:,0]
        return;