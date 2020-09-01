#!/usr/bin/python

"""
version 1.0
author: Felix Plasser
usage: Print out the value of a specified internal coordinate.
"""

import os, sys, locale
import numpy
#import openbabel
import pybel
import file_handler, struc_linalg

def print_info():
    print 'Usage example: int_coor_multi.py dist 4 5 struc.xyz'
    print 'Supported internal coordinates: dist, bend, tors'
    print 'Arguments: -type <xyz, tmol, mol ...>'
    print '   -dig <number of digits in print out>'
    sys.exit()

if len(sys.argv) < 5:
   print_info()

# defaults
file_type = 'xyz'
digits = 4 # output digits

# read command line input
files = []
coors = []  # which coordinates are read out.
args = sys.argv[1:]
while len(args) > 0:
    arg = args.pop(0)

    if arg == 'dist':
        coors += [['dist',int(args.pop(0)), int(args.pop(0))]]
    elif arg == 'bend':
        coors += [['bend',int(args.pop(0)), int(args.pop(0)), int(args.pop(0))]]
    elif arg == 'tors':
        coors += [['tors',int(args.pop(0)), int(args.pop(0)), int(args.pop(0)), int(args.pop(0))]]
    elif arg == '-type':
        file_type = args.pop(0)
    elif arg == '-dig':
        digits = int(args.pop(0))
    elif arg == '-h':
        print_info()
    elif arg[0] == '-':
        print 'Unsupported option: ' + arg
        print 'int_coor.py -h for more information'
        sys.exit()
    else:
        files += [arg]
        

#print coors, files

struc = struc_linalg.structure()
#obconversion = openbabel.OBConversion()
#obconversion.SetInFormat(file_type)
#str_mol = openbabel.OBMol()

tm = file_handler.table_maker(col_widths=[20]+len(coors)*[5+digits])
tm.print_line(['File']+[coor[0] for coor in coors])
for file in files:
  for mol in pybel.readfile(file_type, file):
#  for i in xrange(5):
#    obconversion.ReadFile(str_mol, file)
    struc.get_mol(mol=mol.OBMol,file_path=file, file_type=file_type)
    out_line = [file]
    for coor in coors:
        if coor[0] == 'dist':
            dist = struc.ret_bond_length(*coor[1:])
            out_line += [locale.format("%.*f", (digits, dist))]
        elif coor[0] == 'bend':
            bend = struc.ret_bend(*coor[1:])
            out_line += [locale.format("%.*f", (digits, bend))]
        elif coor[0] == 'tors':
            tors = struc.ret_tors(*coor[1:])
            out_line += [locale.format("%.*f", (digits, tors))]
    tm.print_line(out_line)
#print tm.return_table()
