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

# Optimized parameters
beta = 9.5
u0 = -3.367
T_inf = 2730
chi_len = 4.01
c0 = 0.3
ep = 10
s_y = 0.85

for t, it in enumerate(time):
  TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
  ny = int(TAU[0]) + 1
  nx = int(len(TAU)/ny)
  TAU = TAU.reshape(nx, ny)
  TAU = TAU[1:, 1:]
  
  tau[it] = s_y*np.mean(TAU)


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

# plot final response
fig = plt.figure(figsize=(11, 8.5))
textStr = '\n'.join((
    r'$\beta$ = %.4f eV$^{-1}$' % (beta, ),
    r'$u_0$ = %.4f eV' % (u0, ),
    r'$T_\infty$ = %d K' % (T_inf, ),
    r'$l_\chi$ = %.3f $\AA$' % (chi_len, ),
    r'$c_0$ = %.3f --' % (c0, ),
    r'$\epsilon_0$ = %.3f --' % (ep, ),
    r'$s_y$ = %.3f GPa' % (s_y, )))

ax = plt.subplot(334)
ax.axis('off')
ax.text(0, 0.5, textStr, verticalalignment='center')

# Plot initial CM effective temperature
fileName = fileDir + '/' + 'tem.' + '0'
T0 = np.fromfile(fileName, dtype=np.float32)
ny = int(T0[0]) + 1
nx = int(len(T0)/ny)
T0 = T0.reshape(nx, ny)
T0 = T0[1:, 1:]

plt.subplot(333)
plt.title(r'$T_{0}$')
plt.contourf(T0, cmap='viridis')
plt.colorbar()
plt.axis('off')

lev1 = np.arange(0,2,0.2)
lev2 = np.arange(800,3200,200)

# Plot Final CM effective temperature
fileName = fileDir + '/' + 'tem.' + '{}'.format(100)
Tf_CM = np.fromfile(fileName, dtype=np.float32)
Tf_CM = Tf_CM.reshape(nx, ny)
Tf_CM = Tf_CM[1:, 1:]

plt.subplot(335)
plt.title(r'${T_{f}}^{CM}$')
plt.contourf(Tf_CM, levels=lev2, cmap='viridis')
plt.colorbar()
plt.axis('off')

# Plot Final CM strain field
fileName = fileDir + '/' + 'Exy.' + '{}'.format(100)
Exy_CM = np.fromfile(fileName, dtype=np.float32)
Exy_CM = Exy_CM.reshape(nx, ny)
Exy_CM = 2*Exy_CM[1:, 1:]   # the output of MD is Γ while the output of CM is E, so we need to multiply with 2.

plt.subplot(336)
plt.title(r'${\Gamma_{f}}^{CM}$')
plt.contourf(Exy_CM, levels=lev1, cmap='viridis')
plt.colorbar()
plt.axis('off')

# plot stress-strain curve
plt.subplot(332)
plt.title(r'$\tau(\gamma)$')
plt.plot(strain, tau, 'k--', label='CM')
plt.plot(strain_md, tau_md, 'k:', label='MD')
plt.xlabel(r'$\gamma$ [--]')
plt.ylabel(r'$\tau$ [GPa]')
plt.legend()

# Plot Final MD reference effective Temperature Field
fileName = refDir + 'pe.MD.' + '{}'.format(100)
Uf_MD = np.loadtxt(fileName, skiprows=1)
Tf_Ref = beta*(Uf_MD-u0)*21000
 
plt.subplot(338)
plt.title(r'${T_{f}}^{Ref}$')
plt.contourf(Tf_Ref, levels=lev2, cmap='viridis')
plt.colorbar()
plt.axis('off')
 
# Plot Final MD strain field
fileName = refDir + 'Exy.MD.' + '{}'.format(100)
Exy_Ref = np.loadtxt(fileName, skiprows=1)
 
plt.subplot(339)
plt.title(r'${\Gamma_{f}}^{Ref}$')
plt.contourf(Exy_Ref, levels=lev1, cmap='viridis')
plt.colorbar()
plt.axis('off')

fig.tight_layout()
figName = 'Response_final.png'
plt.savefig(figName)










