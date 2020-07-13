# computeStress.py
# Darius Alix-Williams


import numpy as np
import matplotlib.pyplot as plt
import os, sys, fnmatch

Ez = 1.809642
kB = 8.617330350e-5

beta = 2.22651354472
u0 = -3.36578382182

beta = 10.6832
u0 = -3.366220

# List all files in sct.out directory
fileDir = 'sct_q.out/'
time = np.arange(0, 101)
tau = np.zeros(time.shape)
strain = np.linspace(0, 0.5, len(tau))

tau_md = np.loadtxt("stress_rand2_matlab.txt", delimiter="\t")
tau_md = -1*tau_md[:, 6]/10000.
strain_md = np.linspace(0, 0.5, len(tau_md))

time_md = np.arange(0, len(tau_md))

strain_inc_md = np.linspace(0, .5, len(tau_md))
print(strain_inc_md)
#sys.exit(0)

strain_inc_cm = np.linspace(0, 0.5, len(time))
print(strain_inc_cm)

istrain = np.isin(strain_inc_md, strain_inc_cm)
print(istrain)

#sys.exit(0)
print(np.isin(strain_inc_md, strain_inc_cm))
print(strain_inc_md[istrain])
print(strain_inc_cm)
sys.exit(0)

for t, it in enumerate(time):
    TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
    ny = int(TAU[0]) + 1
    nx = int(len(TAU)/ny)
    TAU = TAU.reshape(nx, ny)
    TAU = TAU[1:, 1:]

    tau[it] = np.mean(TAU)
    

fig = plt.figure()
plt.plot(strain, tau, label='CM--#')
plt.plot(strain_md, tau_md, label='MD')
plt.xlabel('Strain [--]')
plt.ylabel('Stress [Gpa]')
plt.legend()
plt.savefig('stress-strain.png', dpi=200)
plt.close()

sys.exit(0)

# Plot initial, final temperature fields
for timestep in ['0', '100']:
    
    TEM = np.fromfile(fileDir + 'tem.' + timestep, dtype=np.float32)
    ny = int(TEM[0]) + 1
    nx = int(len(TEM)/ny)
    TEM = TEM.reshape(nx, ny)
    TEM = TEM[1:, 1:]
    
    fig = plt.figure()
    plt.contourf(TEM)
    plt.title('CM - ' + timestep)
    plt.colorbar()
    plt.savefig('tem--CM--' + timestep + '.png', dpi=200)
    plt.close()

sys.exit(0)

fileList = os.listdir(fileDir)
pattern = "Exy.*"
#pattern = "X.*"
#pattern = "dev.*"
# pattern = "tem.*"
for pattern in ["tem.*", "Exy.*"]:
    for entry in fileList:
        if fnmatch.fnmatch(entry, pattern):
            # Get file data
            fileData = np.fromfile(fileDir+entry, dtype=np.float32)
            print(fileData)
            print(fileData.shape)
            ny = int(fileData[0]) + 1
            nx = int(len(fileData)/(ny))
            print(nx)
            print(ny)
            #sys.exit()
            fileData = fileData.reshape(nx, ny)
            fileData = fileData[1:, 1:]
            # Plot as a contour
            fig = plt.figure()
            plt.contourf(fileData)
            plt.title(entry)
            plt.colorbar()
            plt.savefig(entry + '.png')
            plt.close()

            print(np.mean(fileData))
