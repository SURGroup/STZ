import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
#plt.switch_backend('Qt5Agg')
#plt.ion()
import os, sys, fnmatch, subprocess

# Parameters and constants
Ez = 1.809642
TZ = 21000
kB = 8.617330350e-5

# Default values - Calibration parameters
beta=8.0
u0=-3.363
 
refDir="/Users/katianakontolati/Projects/continuum_optimization/STZ/CuZr_Ref/"
fileDir = 'sct_q.out/'
PE_0 = refDir+'pe.MD.0'

time = np.arange(0, 101)
tau = np.zeros(time.shape)
strain = np.linspace(0, 0.5, len(tau))
tau_md = np.loadtxt(refDir+"stress_rand2_matlab.txt", delimiter="\t")
tau_md = -1*tau_md[:, 6]/10000.
strain_md = np.linspace(0, 0.5, len(tau_md))
strain_inc_md = np.linspace(0, .5, len(tau_md))
strain_inc_cm = np.linspace(0, 0.5, len(time))
istrain = np.isin(strain_inc_md, strain_inc_cm)
#print(istrain.shape,istrain)
#print((istrain==True).sum())
#print(strain_inc_md[istrain])
#print(tau_md[istrain])
#print('The length of tau_md is',len(tau_md),len(tau))

U0_MD = np.loadtxt(refDir + 'pe.MD.0', skiprows=1)
u_max = np.amax(U0_MD)


# Load final MD effective temperature, compute upper-limiting value
Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
Tf_MD = beta*(Uf_MD-u0)*TZ
T_inf = np.amax(Tf_MD)
u0=u_max-0.0035
u0=u_max-0.004
count = 0

# Reference - Calibration parameters
beta = 9.5
#betas = np.arange(1.0, 13, 1).tolist() # loop
u0 = -3.367
#u0s = np.airange(-3.366, -3.360, 0.001).tolist()
Delta = 2000   # good value
#Deltas = np.arange(1000,10000,600).tolist()
Omega = 473    # good value
#Omegas = np.arange(50,1000,60).tolist() # parameter omega increases the length of the elastic region
chi_len = 36
#chi_lens = np.arange(2,18,2).tolist()
theta = 150   # good value
#thetas = np.arange(50,300,30).tolist()
eps0 = 0.3
#eps0s = np.arange(0.2,4,0.05).tolist()
c0 = 0.3
#c0s = np.arange(0.05,3.5,0.5).tolist()
Omegaeps0 = 100

# Other parameters
mu = 20
#mus = np.arange(10,22,2).tolist()
kappa = 1
sy = 0.85
#syss = np.arange(0.3,1,0.1).tolist()
rho = 7234
#rhos = np.arange(6000,10000,500).tolist()
K = 143.23
#Ks = np.arange(142,143,1).tolist()

# Figure to store stress-strain evolution with kappa variation
f1, ax1 = plt.subplots(1,1)
f2, ax2 = plt.subplots(1,1)
stat = True

Tf_MD = beta*(Uf_MD-u0)*TZ
T_inf = np.amax(Tf_MD)

#command = ['./shear_data_Adam', 'qs', 'null']


for t, it in enumerate(time):
  TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
  ny = int(TAU[0]) + 1
  nx = int(len(TAU)/ny)
  TAU = TAU.reshape(nx, ny)
  TAU = TAU[1:, 1:]
  
  tau[it] = 0.85*np.mean(TAU)


#plt.title(r'${T_{f}}^{CM}$')
#plt.contourf(Tf_CM, cmap='Greys')
#plt.colorbar()
#plt.axis('off')

# Final CM plastic strain
#Exy_CM = np.fromfile(fileDir + 'Exy.100', dtype=np.float32)
#Exy_CM = Exy_CM.reshape(nx, ny)
#Exy_CM = 2.*Exy_CM[1:, 1:]

#ax2 = plt.subplot(336)
#plt.title(r'${\Gamma_{f}}^{CM}$')
#plt.contourf(Exy_CM, cmap='Greys')
#plt.colorbar()
#plt.axis('off')
#f2.savefig('Final_CM_plastic_strain.png')        
"""
# Final MD effective temperature
Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
Tf_MD = beta*(Uf_MD-u0)*TZ
plt.subplot(338)
plt.title(r'${T_{f}}^{MD}$')
plt.contourf(Tf_MD, cmap='Greys')
plt.colorbar()
plt.axis('off')
"""
# Final MD plastic strain
#Exy_MD = np.loadtxt(refDir + 'Exy.MD.100', skiprows=1)
#plt.subplot(339)
#plt.title(r'${\Gamma_{f}}^{MD}$')
#plt.contourf(Exy_MD, cmap='Greys')
#plt.colorbar()
#plt.axis('off')

# Plot Stress-Strain
if stat == True:
    ax1.plot(strain_md, tau_md, color='black', label='MD')
    ax1.set_xlabel('Strain')
    ax1.set_ylabel('Stress')
ax1.plot(strain, tau, color='red', label=beta)
#ax1.set_title('K')
#ax1.legend(title='K',loc='right')
f1.savefig('stress_v_strain.png')    

CHI = np.fromfile(fileDir + 'Exy.' + '100', dtype=np.float32)
ny = int(CHI[0]) + 1
nx = int(len(CHI)/ny)
CHI = CHI.reshape(nx, ny)
CHI = 2*CHI[1:, 1:]
plt.title(r'${\Gamma_{f}}^{CM}$')
plt.contourf(CHI, cmap='viridis')
plt.colorbar()
plt.axis('off')
f2.savefig('chi_field.png')

# Compute strain field distance

def compareMatrices(A, B, verbose=False):
    UA, SA, VA = np.linalg.svd(A, full_matrices=True)
    UB, SB, VB = np.linalg.svd(B, full_matrices=True)

    SA = np.diag(SA)
    SB = np.diag(SB)

    stol=1e-5
    rA = np.linalg.matrix_rank(SA, tol=stol)
    rB = np.linalg.matrix_rank(SB, tol=stol)
    r = np.amax([rA, rB])

    UA = UA[:, :r]
    UB = UB[:, :r]

    R = np.dot(UA.T, UB)
    UR, SR, VR = np.linalg.svd(R, full_matrices=False)
    SR = np.round_(SR, 6)
    SR[np.where(SR >1)] = 1.0
    theta = np.arccos(SR)
    dist = np.sqrt(np.sum(theta**2))
    methodtxt = "Geodesic"

    if verbose:
        text1 = "{} distance of {} determined.\n".format(methodtxt, dist)
        text2 = "Original ranks rA = {} and rB = {}\n".format(rA, rB)
        text3 = "Final rank r = {}.\n".format(r)
        print(text1 + text2 + text3)

    return dist, rA, rB


    



strain =  [0.05, 0.06, 0.08, 0.09, 0.10, 0.11, 0.12, 0.35, 0.5]
gamOBS = np.linspace(0, 0.5, 101)
A = []
B = []
for val in strain:
    gamREF = np.linspace(0, 0.5, 101)
    idREF = np.argmin(np.abs(gamREF-val))
    Aname = refDir+'Exy.MD.'+'{}'.format(idREF)
    A_ = np.loadtxt(Aname, skiprows=1)
    A.append(A_)

    idOBS = np.argmin(np.abs(gamOBS-val))
    Bname = fileDir+'Exy.'+'{}'.format(idOBS)
    B_ = np.fromfile(Bname, dtype=np.float32)
    nyB = int(B_[0]) + 1
    nxB = int(len(B_)/nyB)
    B_ = B_.reshape(nxB, nyB)
    B_ = 2.*B_[1:, 1:]
    B.append(B_)

distance_temp = np.zeros(len(A))
r0_temp = np.zeros(len(A))
rR_temp = np.zeros(len(A))
for k in range(len(A)):
    distance_temp[k], r0_temp[k], rR_temp[k] = compareMatrices(A[k], B[k])

r0 = r0_temp
rR = rR_temp
distance1 = np.mean(distance_temp) 


# Compute stress-strain curve distance

strain_cm = np.linspace(0, 0.5, len(tau))
stress_cm = tau
distance_i = np.zeros(strain_cm.shape)
for i in range(0, 101):
    ix = np.argmin(np.abs(strain_md-strain_cm[i]))
    distance_i[i] = np.abs(tau_md[ix]-stress_cm[i])

distance_i[15:28] *= 5
distance2 = np.sum(distance_i)/np.count_nonzero(distance_i)

distance = (0.25*distance1 + 0.75*distance2)
print('Distance from strain fields is {}, distance from stress-strain curve is {} and total weighted distance is {}'.format(distance1,distance2,distance))





