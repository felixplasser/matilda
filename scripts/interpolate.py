#!/usr/bin/python

"""
version 1.0.0
author: Felix Plasser
usage: example of performing a linear interpolation between two structures with the package struc_linalg.
"""

import os, sys
from struc_manip import struc_linalg

if not len(sys.argv) == 6:
   print 'Five arguments required.'
   print 'Syntax: python interpolate.py <start_struc> <end_struc> <out_dir> <steps_nr> <file_type>'
   sys.exit()
   

start_file = sys.argv[1]
end_file = sys.argv[2]
out_dir = sys.argv[3]
steps_nr = eval(sys.argv[4])
file_type = sys.argv[5]

try:
    os.makedirs(out_dir)
except OSError:
    print 'Output directory could not be created. It either already exists or you have no writing access.'
    
# for calculations an instance of this class with a default molecule has to be defined
mc = struc_linalg.mol_calc(def_file_path=start_file, file_type=file_type)

# load the structures
st = struc_linalg.structure(name='st') # optional name for output
st.read_file(file_path=start_file, file_type=file_type)

en_temp = struc_linalg.structure(name='en_t')
en_temp.read_file(file_path=end_file, file_type=file_type)

# superimpose the structures
en = en_temp.ret_superimposed_structure(struc=st, mass_wt_pw=1, name='en')

# define the vector for the linear interpolation
diff_vect = mc.subtract(en, st)

# create <steps_nr> structures and put them into separate folders
for i in xrange(1,steps_nr+1):
    x = i * 1. / (steps_nr + 1)
    plus_struc = mc.scalar_mult(x, diff_vect)
    grid_struc = mc.add(st, plus_struc)

    try:
        os.makedirs(out_dir + '/struc_%6.4f'%(x))
    except OSError:
        print 'Output directory could not be created. It either already exists or you have no writing access.'

    grid_struc.make_coord_file(os.path.join(out_dir,  'struc_%6.4f'%(x), 'coord.'+file_type), file_type=file_type) 
