#!/usr/bin/python

"""
version 1.0.0
author: Felix Plasser
usage: Print distances between the initial geometry and every consecutive geometry.
"""

import os, sys, locale
import numpy
import openbabel
import struc_linalg_copy as struc_linalg # read_file_3xN_matrix does not work in the original struc_linalg
  # with this the original object is also changed during superposition

if len(sys.argv) < 2:
   print 'At least one argument required.'
   print 'Syntax: python geom_change.py <file> [-type <type>] [-mwp <mass_weight_power>] [-fit <fitting>] [-dig <digits>] [-ref <ref_struc>] [plot]'
   sys.exit()
   
# defaults
file_type = 'xyz'
mass_wt_pw = 1 # mass weighting to which power (1 is regular mass weighting), only valid if fitting=1
fitting = True # if structures are superimposed
digits = 4 # output digits
ref_file = ''
plot = False

args = sys.argv[1:]
while len(args) > 0:
    arg = args.pop(0)

    if arg == '-type':
        file_type = args.pop(0)
    elif arg == '-mwp':
        mass_wt_pw = eval(args.pop(0))
    elif arg == '-fit':
        fitting = eval(args.pop(0))
    elif arg == '-dig':
        digits = int(args.pop(0))
    elif arg == '-ref':
        ref_file = args.pop(0)
    elif arg == '-plot':
        try:
            import pylab
            plot = True
        except:
            print 'Plotting not possible'
    elif arg[0] == '-':
        print 'Unsupported option: ' + arg
        sys.exit()
    else:
        file = arg


if ref_file == '': ref_file = file
struc0 = struc_linalg.structure()
struc0.read_file(file_path=ref_file, file_type=file_type)

        ## read in the structures with openbabel
distances = []
temp_struc = struc_linalg.structure()
mc = struc_linalg.mol_calc(def_file_path=ref_file, file_type=file_type)

obconversion = openbabel.OBConversion()
obconversion.SetInFormat(file_type)
mol = openbabel.OBMol()

#print 'timestep:'
# the first structure is read in
notatend = obconversion.ReadFile(mol, file)

# the other structures are read in
while notatend:
    temp_struc.get_mol(mol, file_type) # read in the data
    
    # compute the distance to the first structure
    if fitting:
        comp_struc = temp_struc.ret_superimposed_structure(struc0, mass_wt_pw=mass_wt_pw)
        #print comp_struc.ret_vector()
    else:
        comp_struc = temp_struc

    distances.append(mc.distance(struc0, comp_struc, mass_wt_pw=mass_wt_pw))
    #distances.append(mc.distance(struc0, temp_struc, mass_wt_pw=mass_wt_pw))

    mol = openbabel.OBMol()
    notatend = obconversion.Read(mol)

    ## print out results

print 'Distances between refererence structure and consecutive structures'
for i,distance in enumerate(distances):
    print i, locale.format("%.*f", (digits, distance))

if plot:
    pylab.plot(range(len(distances)), distances)
    pylab.show()
