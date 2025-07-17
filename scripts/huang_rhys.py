#!/usr/bin/python3

from __future__ import print_function

"""
version 1.0.0
author: Felix Plasser
usage: Compute Huang-Rhys parameters and related quantities.
"""

import os, sys
import numpy # tmp
from matilda import vib_molden
from matilda import struc_linalg
from matilda import units

print("huang_rhys.py <gs_struc> <es_struc> <vib_mld>")

if len(sys.argv) < 4:
   print('At least 3 arguments required.')
   sys.exit()
   
# read in data.
gs_file = sys.argv[1]
es_file = sys.argv[2]
vib_file = sys.argv[3]

file_type = 'xyz'

# TODO: One could also read gs_struc from vib_mld
gs_struc = struc_linalg.structure()
gs_struc.read_file(file_path=gs_file, file_type=file_type)

es_struc = struc_linalg.structure()
es_struc.read_file(file_path=es_file, file_type=file_type)
# TODO: Could do superposition here

diff_vec = es_struc.ret_vector() - gs_struc.ret_vector()
print("Distance: %.4f Ang"%numpy.sum(diff_vec*diff_vec)**.5)

M = gs_struc.ret_mass_vector(power=0.5, rep=3)
mdiff_vec = diff_vec * M
print("MW Dist: %.5f amu**.5 Ang"%(numpy.sum(mdiff_vec*mdiff_vec))**.5)
print("MW MDist   : %.5f Ang"%(numpy.sum(mdiff_vec*mdiff_vec)/numpy.sum(M*M))**.5)

mc = struc_linalg.mol_calc(gs_file, file_type=file_type)
print("RMSD /3**.5: %.5f Ang"%(mc.RMSD(gs_struc, es_struc) / 3**.5))

# Read information from vibration file
vmol = vib_molden.vib_molden()
vmol.read_molden_file(vib_file)
Kmat = vmol.ret_vib_matrix()

print("Kmat")
print(Kmat)
print(numpy.dot(Kmat, Kmat.T))

Om = vmol.ret_freqs() / units.energy['rcm']

#U_TO_AMU = 1./5.4857990943e-4

#M1 = gs_struc.ret_mass_vector(power=1, rep=3) / U_TO_AMU
#Mm = gs_struc.ret_mass_vector(power=0.5, rep=3) * U_TO_AMU**.5

for mode in Kmat:
    print('mode', mode)
    mode /= M
    norm = numpy.sum(mode * mode)
    if norm > 1E-6:
        mode *= norm**(-.5)
    print('mode2', mode)

#print("KT K")
#%KK = numpy.dot(Kmat, Kmat.T)
#print(KK[-8:,-8:])
#print(numpy.sum(KK*KK))

#MM = mc.ret_mass_matrix(0)
#K2 = numpy.dot(MM, numpy.dot(Kmat, MM))

#print("K2")
#print(numpy.dot(K2.T, K2)[-8:,-8:])

#print("dot")
#KK = numpy.dot(Kmat.T, numpy.dot(MM, Kmat))
#print(KK[-8:,-8:])

#print(Kmat[1])

dQ = Om**.5 * numpy.dot(mdiff_vec, Kmat.T) / units.length['A']

print("mode disp    S_i")
for i, val in enumerate(dQ):
    print("%3i % 7.4f % 7.4f"%(i+1, val, val*val/2.))
