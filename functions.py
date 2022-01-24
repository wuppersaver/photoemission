#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 10:50:04 2022

@author: brunocamino
"""

def calc_surface_energy(bulk_energy,slab_energy,n_units):
    #bulk_energy (float): energy of the bulk
    #slab_energy (float): energy of the slab
    #n_units (int): number of bulk units in the slab
    
    e_surf = 0.5*slab_energy-(bulk_energy*n_units)
    
    return e_surf


def get_qe(optados_output):
    #This is for you to develop :)
    pass