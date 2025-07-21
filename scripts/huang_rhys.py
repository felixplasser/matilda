#!/usr/bin/python3

from __future__ import print_function

"""
version 1.0.0
author: Felix Plasser
usage: Compute Huang-Rhys parameters and related quantities.
"""

import os, sys
import numpy # tmp
import matplotlib.pyplot as plt
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
mdiff_vec *= units.mass['amu']**(-.5) / units.length['A'] # Conversion to atomic units

mc = struc_linalg.mol_calc(gs_file, file_type=file_type)
print("RMSD /3**.5: %.5f Ang"%(mc.RMSD(gs_struc, es_struc) / 3**.5))

# Read information from vibration file
vmol = vib_molden.vib_molden()
vmol.read_molden_file(vib_file)
Kmat = vmol.ret_vib_matrix()

freqs_cm = vmol.ret_freqs()   # Frequencies in cm^-1
Om = freqs_cm / units.energy['rcm']

# TODO: put this into vib_molden
for mode in Kmat:
    mode *= M
    norm = numpy.sum(mode * mode)
    if norm > 1E-8:
        mode *= norm**(-.5)

#print("Final Kmat:")
#print(Kmat.T) -> This seems to work

dQ = Om**.5 * numpy.dot(mdiff_vec, Kmat.T)
S_factors = 0.5 * dQ**2


valid_indices = freqs_cm > 1E-8
freqs_cm = freqs_cm[valid_indices]
Om = Om[valid_indices]
dQ = dQ[valid_indices]
S_factors = S_factors[valid_indices]

# Reorganization energy (Hartree and eV)
CM_TO_AU = 1.0 / 219474.6313705
HARTREE_TO_EV = 27.2114

freqs_au = freqs_cm * CM_TO_AU
reorg_energy_au = numpy.sum(S_factors * freqs_au)
reorg_energy_ev = reorg_energy_au * HARTREE_TO_EV

print("\nMode   Frequency(cm^-1)   dQ     Huang-Rhys S")
for i in range(len(freqs_cm)):
    print(f"{i+1:3d}    {freqs_cm[i]:10.2f}   {dQ[i]:6.4f}   {S_factors[i]:10.6f}")

print(f"\nReorganization Energy: {reorg_energy_ev:.6f} eV")


# Print top 10 modes
top_indices = numpy.argsort(S_factors)[-10:][::-1]
print("\nTop 10 modes with highest Huang-Rhys factors:")
print(f"{'Mode':>4} {'Freq (cm^-1)':>15} {'S_i':>10}")
for idx in top_indices:
    print(f"{idx+1:4d} {freqs_cm[idx]:15.2f} {S_factors[idx]:10.6f}")

# Plot Huang-Rhys Spectrum
plt.figure(figsize=(8,5))
plt.bar(freqs_cm, S_factors, width=5, align='center', color='mediumblue', edgecolor='black', linewidth=0.7)
plt.xlabel('Vibrational frequency (cm$^{-1}$)')
plt.ylabel('Huang-Rhys factor $S_i$')
plt.title('Huang-Rhys Spectrum')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("huang_rhys_spectrum.png", dpi=300)
#plt.show()

