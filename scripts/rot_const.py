#!/usr/bin/env python3

import numpy as np 
import argparse
from matilda.atominfo import Z_exact_mass_dict, symbol_Z_dict

all_element_mass = {element_symbol: Z_exact_mass_dict[element_mass] for element_symbol,element_mass in symbol_Z_dict.items() if element_mass in Z_exact_mass_dict}

parser = argparse.ArgumentParser(prog = 'rot_const',
                                description = 'Computes rotational constants (in MHz) of the molecule from geomery')
parser.add_argument('-f', '--filename' , default = 'orca.xyz')
args = parser.parse_args()

filename = args.filename

def file_read(filename):
    elements = []
    coordinates = []

    with open(filename,'r') as file:
        lines = file.readlines()
        number_of_atoms = int(lines[0])
        for line in lines[2:2+number_of_atoms]:
            words = line.split()
            elements.append(words[0])
            coordinates.append([float(coor) for coor in words[1:4]])

    return elements , np.array(coordinates)

def mass_center(elements,coordinates):
    molecular_mass = sum(all_element_mass[e] for e in elements)
    mass_center_coordinates = sum(all_element_mass[e] * coordinates[c] for c,e in enumerate(elements))
    
    return mass_center_coordinates/molecular_mass

def moment_of_inertia(elements, coordinates):
    coordinates = coordinates * 1.0e-10
    masses = np.array([all_element_mass[e] * 1.66054e-27 for e in elements])

    mc = mass_center(elements, coordinates)
    coordinates -= mc

    I = np.zeros((3,3))
    
    for i in range(len(elements)):
        x ,y ,z = coordinates[i]
        m = masses[i]
        I += m* np.array([
            [y**2+z**2, -x*y, -x*z],
            [-x*y, x**2+z**2, -y*z],
            [-x*z, -y*z, x**2+y**2]
        ])
    
    return I
        

def rot_constant(filename):
    elements, coordinates = file_read(filename)
    I = moment_of_inertia(elements, coordinates)
    eigenvalues,_ = np.linalg.eigh(I)

    h = 6.626e-34
    rot_const_MHz = h/(8 * np.pi**2 * eigenvalues * 1e6)
    print("Rotational constants in MHz :", " ".join(f"{i:12.6f}" for i in rot_const_MHz))
rot_constant(filename)
