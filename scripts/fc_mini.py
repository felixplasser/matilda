#!/usr/bin/python3

from __future__ import print_function

"""
version 1.0.0
author: Felix Plasser
usage: Spectrum using Franck-Condon progression formula based on Huang-Rhys factors.
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

print("Select mode: 1 = absorption (abs) or 2 = emission (emi)")
choice = input("Enter 1 or 2: ").strip()

if choice == "1":
    mode_type = "abs"
elif choice == "2":
    mode_type = "emi"
else:
    print("Invalid input. Please enter 1 or 2.")
    sys.exit()


hr_data_file = sys.argv[1]

try:
    delta_E_hartree = float(input("Electronic adiabatic energy (in hartree): "))
except ValueError:
    print("Invalid input. Please enter a valid floating-point number for energy.")
    sys.exit()

E_0_cm = delta_E_hartree * units.energy['rcm']
print(f"Electronic adiabatic energy: {E_0_cm:.4f} cm^-1")

f = float(input("\nOscillator strength (dimensionless): "))

# input w_min
try:
    w_min = float(input("\nw_min(in cm^-1): "))
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
sigma = sum(0.5 * (omega**2) * S for S, omega in sig_modes) ** 0.5

print()
print(f"Using modes with omega > {w_min: .4f} cm^-1 and S > {S_min: .4f} as FC-active modes.")
print(f"Number of active modes : {len(FC_modes)}")
print(f"Number of sigma modes  : {len(sig_modes)}")
print(f"Computed gaussian sigma: {sigma:.4f} cm^-1")


# Spectrum generation parameters
k_max = int(max(S + 5 * math.sqrt(S) for S, _ in FC_modes))
#print(f"k_max= {k_max}")
n_points = 2000          # Number of points in the energy grid

# Define energy grid (in cm^-1)
fc_width = sum(S * omega for S, omega in FC_modes)

if mode_type == "emi":
    E_min = E_0_cm - 10000
    E_max = E_0_cm + 10000
    #E_min = E_0_cm - fc_width - 20*sigma
    #E_max = E_0_cm + 10*sigma
else:
    E_min = E_0_cm - 5*sigma
    E_max = E_0_cm + fc_width + 50*sigma

#E_min = E_0_cm - 10000
#E_max = E_0_cm + 10000

energy_grid = numpy.linspace(E_min, E_max, n_points)

# Gaussian line shape in cm^-1
def gaussian(x, x0, sigma):
    return numpy.exp(-0.5 * ((x - x0)/sigma)**2) / (sigma * numpy.sqrt(2 * numpy.pi))


def build_spectrum(modes, k_max, E_0_cm, energy_grid, sigma):
    spectrum = numpy.zeros_like(energy_grid)
    S_total = sum(S for S, _ in modes)
    prefactor = math.exp(-S_total)

    stick_energies = []
    stick_intensities = []

    # Recursive loop
    def loop_over_modes(idx, E_shift, intensity):
        if idx == len(modes):
            # All modes handled -> place peak
            E_transition = E_0_cm + E_shift
            spectrum[:] += intensity * gaussian(energy_grid, E_transition, sigma)

            stick_energies.append(E_transition)
            stick_intensities.append(intensity)
            return
        S, omega = modes[idx]
        for k in range(k_max + 1):
            new_intensity = intensity * (S**k) / math.factorial(k)
            if mode_type == "emi":
                new_Eshift = E_shift - k * omega
            else:
                new_Eshift = E_shift + k * omega
            loop_over_modes(idx + 1, new_Eshift, new_intensity)

    # Start recursion
    loop_over_modes(0, 0.0, prefactor)
    return spectrum, stick_energies, stick_intensities


spectrum, stick_energies, stick_intensities = build_spectrum(FC_modes, k_max, E_0_cm, energy_grid, sigma)

# Emission
if mode_type == "emi":
    ein_coeff_emi = (10000 * 4 * math.pi * units.constants['h']) / (units.mass['kg'] * units.constants['c0']) * f * (E_0_cm**2)
    Lamda = 1/ein_coeff_emi
    print("\nLineshapes are ignored")
    print(f"Einstein coefficient of Spontaneous Emission A: {ein_coeff_emi:.6e} s^-1")
    print(f"Excited state lifetime:{Lamda:.6e} s")

    spectrum /= spectrum.max()

    stick_intensities_normalised = numpy.array(stick_intensities) / max(stick_intensities)

    wavelength_grid = 1e7 / energy_grid
    electronvolts_grid = energy_grid / 8065.54

    stick_wavelengths = 1e7 / numpy.array(stick_energies)
    stick_ev = numpy.array(stick_energies) /8065.54

    # Save emission data
    output_file = "vibronic_emission_data.txt"
    with open(output_file, "w") as f:
        f.write("Wavenumber(cm^-1)   Wavelength(nm)   Electronvolts(eV)   Intensity(normalized)\n")
        for E, L, V, I in zip(energy_grid, wavelength_grid, electronvolts_grid, spectrum):
            f.write(f"{E:15.6f}  {L:15.6f}  {V:15.6f}  {I:15.8f}\n")

    print()
    print(f"Emission spectrum saved to '{output_file}'")

    # Save stick spectrum data
    stick_output_file = "vibronic_emission_stick_data.txt"
    with open(stick_output_file, "w") as f:
        f.write("Wavenumber(cm^-1)   Wavelength(nm)   Electronvolts(eV)   Intensity(normalized)\n")
        for E, L, V, I in zip(stick_energies, stick_wavelengths, stick_ev, stick_intensities_normalised):
            f.write(f"{E:15.6f}  {L:15.6f}  {V:15.6f}  {I:15.8f}\n")
    print(f"Stick spectrum saved to '{stick_output_file}'")


    # Plot emission
    plt.figure(figsize=(8,5))
    plt.plot(energy_grid, spectrum, color='darkred')
    plt.vlines(stick_energies, 0, stick_intensities_normalised, color='blue', linewidth=1, label='Stick spectrum')
    plt.xlim(E_min, E_max)
    plt.xlabel('Wavenumber (cm$^{-1}$)')
    plt.ylabel('Normalized Intensity')
    plt.title('Simulated Vibronic Emission Spectrum')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig("vibronic_emission.png", dpi=300)
    print("Plot saved as 'vibronic_emission.png'")

    # Wavelength plot
    plt.figure(figsize=(8,5))
    plt.plot(wavelength_grid, spectrum, color='orange')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Normalized Intensity')
    plt.title('Simulated Vibronic Emission Spectrum (Wavelength)')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig("vibronic_emission_nm.png", dpi=300)
    print("Plot saved as 'vibronic_emission_nm.png'")

    #Electronvolts plot
    plt.figure(figsize=(8,5))
    plt.plot(electronvolts_grid, spectrum, color='blue')
    plt.xlabel('Electronvolts (eV)')
    plt.ylabel('Normalized Intensity')
    plt.title('Simulated Vibronic Emission Spectrum (Electronvolts)')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig("vibronic_emission_ev.png", dpi=300)
    print("Plot saved as 'vibronic_emission_ev.png'")

    sys.exit()

# Absorption
ein_coeff_abs = math.pi / (1000 * units.mass['kg'] * 100 * units.constants['c'] * units.constants['c0']) * f / E_0_cm
print("\nLineshapes are ignored and, in the electric dipole approximation we can obtain Einstein coefficient of absorption(B).")
print(f"B: {ein_coeff_abs:.6e} s gm^-1")

print("\nSelect absorption cross-section treatment:")
print("1 = Assume omega/omega_I0 = 1 sharp lineshape approximation")
print("2 = Use full omega/omega_I0 frequency dependence")

pref_choice = input("Enter 1 or 2: ").strip()

if pref_choice not in ["1", "2"]:
    print("Invalid choice. Please restart and enter 1 or 2.")
    sys.exit()

cross_sec_constant = (100 * units.constants['h'])/(units.mass['kg'] * units.constants['c'] * (8 * math.pi)**0.5 * units.constants['c0'])
epsilon_prefactor = cross_sec_constant * unit.constants['Nl'] / math.log(10)

print()

if pref_choice == "1":
    print("Using Option 1:")
    print("Assuming omega/omega_I0 = 1 (sharp lineshape approximation).")
    print("The absorption cross section has uniform scaling across the spectrum.")

    cross_sec_sigma = cross_sec_constant * f / sigma
    print(f"\nCharacteristic peak absorption cross section: {cross_sec_sigma:.6e} cm^2")

    epsilon_spectrum = epsilon_prefactor * spectrum

elif pref_choice == "2":
    print("Using Option 2:")
    print("Including full frequency-dependent factor (omega/omega_I0).")
    print("The absorption cross section now varies at each frequency point.")
    print("No single fixed peak formula applies.")

    omega_ratio = energy_grid / E_0_cm
    epsilon_spectrum = epsilon_prefactor * omega_ratio * spectrum


wavelength_grid = 1e7 / energy_grid
electronvolts_grid = energy_grid / 8065.54

output_file = "vibronic_spectrum_data.txt"
with open(output_file, "w") as f:
    f.write("Wavenumber(cm^-1)   Wavelength(nm)   Electronvolts(eV)   Molar Extinction Coefficients(M^-1 cm^-1)\n")
    for E, L, V, I in zip(energy_grid, wavelength_grid, electronvolts_grid, epsilon_spectrum):
        f.write(f"{E:15.6f}  {L:15.6f}  {V:15.6f}  {I:15.8f}\n")

print()
print(f"Absorption spectrum saved to '{output_file}'")

# Plot absorption
plt.figure(figsize=(8,5))
plt.plot(energy_grid, epsilon_spectrum, color='darkgreen')
plt.xlabel('Wavenumber (cm$^{-1}$)')
plt.ylabel('Molar Extinction Coefficient (M$^{-1}$ cm$^{-1}$)')
plt.title('Simulated Vibronic Absorption Spectrum')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("vibronic_spectrum.png", dpi=300)
print("Plot saved as 'vibronic_spectrum.png'")

# Wavelength plot
plt.figure(figsize=(8,5))
plt.plot(wavelength_grid, epsilon_spectrum, color='purple')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Molar Extinction Coefficient (M$^{-1}$ cm$^{-1}$)')
plt.title('Simulated Vibronic Absorption Spectrum (Wavelength)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("vibronic_spectrum_nm.png", dpi=300)
print("Plot saved as 'vibronic_spectrum_nm.png'")

# Electronvolts plot
plt.figure(figsize=(8,5))
plt.plot(electronvolts_grid, epsilon_spectrum, color='brown')
plt.xlabel('Electronvolts (eV)')
plt.ylabel('Molar Extinction Coefficient (M$^{-1}$ cm$^{-1}$)')
plt.title('Simulated Vibronic Absorption Spectrum (Electronvolts)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("vibronic_spectrum_ev.png", dpi=300)
print("Plot saved as 'vibronic_spectrum_ev.png'")