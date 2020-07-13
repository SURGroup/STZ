# computeStress.py
# Darius Alix-Williams (d.alixwill@gmail.com)

import numpy as np
import matplotlib.pyplot as plt
import os, sys, fnmatch

# Constants
min_strain = 0
max_strain = 0.5
nSnapshots = 101

# List all files in sct.out directory
fileDir = 'sct_q.out/'
refOutFileName = "stress_rand2_matlab.txt"

# Get reference stress, strain from MD sourcefile
stress_md = np.loadtxt(refOutFileName, delimiter="\t")
stress_md = -1*stress_md[:, 6]/10000.
nPoints = len(stress_md)
strain_md = np.linspace(min_strain, max_strain, len(stress_md))

# Get continuum stress, strain from output snapshots
time = np.arange(0, nSnapshots)
stress_cm = np.zeros(time.shape)
strain_cm = np.linspace(min_strain, max_strain, nSnapshots)
istrain = np.isin(strain_md, strain_cm)

for t, it in enumerate(time):
    TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
    ny = int(TAU[0]) + 1
    nx = int(len(TAU)/ny)
    TAU = TAU.reshape(nx, ny)
    TAU = TAU[1:, 1:]
    stress_cm[it] = np.mean(TAU)
    
# Plot the result
fig = plt.figure()
plt.plot(strain_cm, stress_cm, label='CM')
plt.plot(strain_md, stress_md, label='MD')
plt.xlabel('Strain [--]')
plt.ylabel('Stress [Gpa]')
plt.legend()
plt.savefig('stress-strain.png', dpi=200)
plt.close()
