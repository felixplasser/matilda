#!/usr/bin/env python3

"""
version 1.0
author: Felix Plasser
usage: Print out the value of a specified internal coordinate.
"""
import argparse
import numpy, sys
from openbabel import pybel
from matilda import file_handler, struc_linalg

parser = argparse.ArgumentParser(
    prog='int_coor_multi.py',
    description='Outputs the value of the speciefied internal coordinate'
)
parser.add_argument('-d', '--dist', nargs = 2, type = int, default = [1,2], help = 'The indices of atoms for distance')
parser.add_argument('-b', '--bend', nargs = 3, type = int, default = [1 ,2 ,4], help = 'The indices of atoms for bend angle')
parser.add_argument('-t', '--tors', nargs = 4, type = int, default = [5, 1 ,2 ,6], help = 'The indices of atoms for dihedral angle')
parser.add_argument('-f', '--filename', default = ['dyn.xyz'])
parser.add_argument('-ft', '--filetype', default = 'xyz')
parser.add_argument('-dg', '--digits',type = int,  default = 4, help = 'Number of decimal points')
args = parser.parse_args()

dist = args.dist
bend = args.bend
tors = args.tors
files = args.filename
file_type = args.filetype
digits = args.digits # output digits

def print_info():
    print('Usage example: int_coor_multi.py dist 4 5 struc.xyz')
    print('Supported internal coordinates: dist, bend, tors')
    print('Arguments: -type <xyz, tmol, mol ...>')
    print('   -dig <number of digits in print out>')
    sys.exit()

coors = []  # which coordinates are read out.
if args.dist:
  coors.append(['dist'] + args.dist)
if args.bend:
    coors.append(['bend'] + args.bend)
if args.tors:
    coors.append(['tors'] + args.tors)

struc = struc_linalg.structure()

tm = file_handler.table_maker(col_widths=[20]+len(coors)*[5+digits])
tm.print_line(['File']+[coor[0] for coor in coors])
for file in files:
  for mol in pybel.readfile(file_type, file):
    struc.get_mol(mol=mol.OBMol,file_path=file, file_type=file_type)
    out_line = [file]
    for coor in coors:
        if coor[0] == 'dist':
            dist = struc.ret_bond_length(*coor[1:])
            out_line += ["%.*f"%(digits, dist)]
        elif coor[0] == 'bend':
            bend = struc.ret_bend(*coor[1:])
            out_line += ["%.*f"%(digits,bend)]
        elif coor[0] == 'tors':
            tors = struc.ret_tors(*coor[1:])
            out_line += ["%.*f"%(digits,tors)]
    tm.print_line(out_line)
