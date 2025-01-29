# This script plots the Mach and pitch scans for each of the relevant aerodynamic parameters

import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science'])
plt.style.use(['no-latex'])
# Read in the .dat files, excluding the first column


mach_data = np.loadtxt("./mach_scan/all_aero_coefs.dat", usecols=range(1, 22))

pitch_data = np.loadtxt("./aoa_scan/all_aero_coefs.dat", usecols=range(1, 22))

# First plot lift to drag ratio vs angle of attack
plt.figure(figsize=(5,5))
ax = plt.gca()

params = {
    'axes.labelsize': 8,
    'font.size': 8,
    'font.family': 'serif',
    'legend.fontsize': 8,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    #'text.usetex': True,
}

plt.rcParams.update(params)

# set color palette
# colors = plt.cm.viridis(np.linspace(0, 1, 7))

plt.plot(pitch_data[:, 0], pitch_data[:, 8], label=r'$C_L/C_D$')

# title and labels
plt.title('Lift to Drag Ratio vs Angle of Attack')
plt.xlabel('Angle of Attack (deg)')
plt.ylabel(r'$C_L/C_D$')

plt.savefig("./aoa_scan/lift_to_drag.pdf")

# Plot lift, drag, and moment coefficients vs angle of attack
plt.figure(figsize=(5,5))
ax = plt.gca()
params = {
    'axes.labelsize': 8,
    'font.size': 8,
    'font.family': 'serif',
    'legend.fontsize': 8,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    #'text.usetex': True,
}
plt.rcParams.update(params)
# set color palette
# colors = plt.cm.viridis(np.linspace(0, 1, 7))
plt.plot(pitch_data[:, 0], pitch_data[:, 5], label=r'$C_D$')
plt.plot(pitch_data[:, 0], pitch_data[:, 7], label=r'$C_L$')
plt.plot(pitch_data[:, 0], pitch_data[:, 13], label=r'$C_{M_{\alpha}}$')

# title and labels
plt.title('Aerodynamic Coefficients vs Angle of Attack')
plt.xlabel('Angle of Attack (deg)')
plt.ylabel('Coefficient')
plt.legend(loc='upper left')
plt.savefig("./aoa_scan/aero_coefs.pdf")

# Plot drag coefficient vs mach number
plt.figure(figsize=(5,5))
ax = plt.gca()
params = {
    'axes.labelsize': 8,
    'font.size': 8,
    'font.family': 'serif',
    'legend.fontsize': 8,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    #'text.usetex': True,
}
plt.rcParams.update(params)
# set color palette
# colors = plt.cm.viridis(np.linspace(0, 1, 7))
plt.plot(mach_data[:, 0], mach_data[:, 5], label=r'$C_D$')
# title and labels
plt.title('Drag Coefficient vs Mach Number')
plt.xlabel('Mach Number')
plt.ylabel(r'$C_D$')
plt.savefig("./mach_scan/drag.pdf")



