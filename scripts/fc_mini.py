#!/usr/bin/python3

from __future__ import print_function

"""
version 1.0.0
author: Felix Plasser
usage: Spectrum using Franck-Condon progression formula based on Huan-Rhys factors.
"""

import os, sys
import math
import numpy
import matplotlib.pyplot as plt
from matilda import units

print("fc_mini <hr_data>")

if len(sys.argv) < 2:
    print("Error: hr_data.txt file required")
    sys.exit()

hr_data_file = sys.argv[1]

try:
    delta_E_hartree = float(input("Electronic adiabatic energy (in hartree): "))
except ValueError:
    print("Invalid input. Please enter a valid floating-point number for energy.")
    sys.exit()

E_0_cm = delta_E_hartree * units.energy['rcm']
print(f"\nElectronic adiabatic energy: {E_0_cm:.4f} cm^-1")


# input w_min
try:
    w_min = float(input("w_min(in cm^-1): "))
except ValueError:
    print("Invalid input. Please enter a number.")
    sys.exit()

# input S_min
try:
    S_min = float(input("          S_min: "))
except ValueError:
    print("Invalid input. Please enter a number.")
    sys.exit()


# Extract freqs and S factors
freqs_cm = []
S_factors = []

for line in open(hr_data_file, 'r').readlines()[1:]:  # Skip header
    parts = line.strip().split()
    if len(parts) < 3:
        continue
    freq = float(parts[1])
    S = float(parts[2])
    freqs_cm.append(freq)
    S_factors.append(S)

if not freqs_cm or not S_factors:
    print("Error: No valid data found in hr_data.txt")
    sys.exit()

# Relevant modes
sig_modes = []
FC_modes = []

for imode in range(len(freqs_cm)):
    omega = freqs_cm[imode]
    S     = S_factors[imode]
    if S >= S_min and omega >= w_min:
        FC_modes += [(S, omega)]
    else:
        sig_modes += [(S, omega)]

if len(FC_modes) == 0:
    print("No modes above the given S_min and w_min for FC progression")
    sys.exit()

# --- compute sigma from Eq. (8) sigma = sqrt(0.5 * sum_i (omega_i**2 * S_i)) ---
sigma = 0.5 * sum((omega**2) * S for S, omega in sig_modes) ** 0.5

print()
print(f"Using modes with omega > {w_min: .4f} cm^-1 and S > {S_min: .4f} as FC-active modes.")
print(f"Number of active modes : {len(FC_modes)}")
print(f"Number of sigma modes  : {len(sig_modes)}")
print(f"Computed gaussian sigma: {sigma:.4f} cm^-1")


# Spectrum generation parameters
k_max = 5               # Max vibrational quantum number for progression
n_points = 2000         # Number of points in the energy grid


# Define energy grid (in cm^-1)
E_min = E_0_cm - 1000
E_max = E_0_cm + 10000
energy_grid = numpy.linspace(E_min, E_max, n_points)


# Gaussian line shape in cm^-1
def gaussian(x, x0, sigma):
    return numpy.exp(-0.5 * ((x - x0)/sigma)**2) / (sigma * numpy.sqrt(2 * numpy.pi))


def build_spectrum(modes, k_max, E_0_cm, energy_grid, sigma):
    spectrum = numpy.zeros_like(energy_grid)
    S_total = sum(S for S, _ in modes)
    prefactor = math.exp(-S_total)

    # Recursive loop
    def loop_over_modes(idx, E_shift, intensity):
        if idx == len(modes):
            # All modes handled -> place peak
            E_transition = E_0_cm + E_shift
            spectrum[:] += intensity * gaussian(energy_grid, E_transition, sigma)
            return
        S, omega = modes[idx]
        for k in range(k_max + 1):
            new_intensity = intensity * (S**k) / math.factorial(k)
            new_Eshift = E_shift + k * omega
            loop_over_modes(idx + 1, new_Eshift, new_intensity)

    # Start recursion
    loop_over_modes(0, 0.0, prefactor)
    return spectrum

spectrum = build_spectrum(FC_modes, k_max, E_0_cm, energy_grid, sigma)

spectrum /= spectrum.max()

# Plot the vibronic spectrum in cm^-1
plt.figure(figsize=(8,5))
plt.plot(energy_grid, spectrum, color='darkgreen')
plt.xlabel('Wavenumber (cm$^{-1}$)')
plt.ylabel('Normalized Intensity')
plt.title('Simulated Vibronic Spectrum')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("vibronic_spectrum.png", dpi=300)
print("Spectrum saved to 'vibronic_spectrum.png'")

