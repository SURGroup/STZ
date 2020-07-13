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


beta=8.0
u0=-3.3626
 
theta = 100
sy = 0.85
mu = 20
K = 143.23
rho = 7234
c0 = 0.3
"""
l_chi -> [1.0, 5.0], default 2 Angstroms
kappa -> [0.009, 0.02], default ~ 0.012 --
"""

refDir="/home-1/dwill164@jhu.edu/projects/04--Parameter_Optimization/CuZr_Ref/"
refDir="/home-3/kkontol1@jhu.edu/CuZr_Ref/"
fileDir = 'sct_q.out/'
PE_0 = refDir+'pe.MD.0'

time = np.arange(0, 101)
tau = np.zeros(time.shape)
q = np.zeros(time.shape)
s = np.zeros(time.shape)
dev = np.zeros(time.shape)
strain = np.linspace(0, 0.5, len(tau))
tau_md = np.loadtxt(refDir+"stress_rand2_matlab.txt", delimiter="\t")
tau_md = -1*tau_md[:, 6]/10000.
strain_md = np.linspace(0, 0.5, len(tau_md))
strain_inc_md = np.linspace(0, .5, len(tau_md))
strain_inc_cm = np.linspace(0, 0.5, len(time))
istrain = np.isin(strain_inc_md, strain_inc_cm)

U0_MD = np.loadtxt(refDir + 'pe.MD.0', skiprows=1)
u_max = np.amax(U0_MD)


# Load final MD effective temperature, compute upper-limiting value
Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
Tf_MD = beta*(Uf_MD-u0)*TZ
T_inf = np.amax(Tf_MD)
u0=u_max-0.0035
u0=u_max-0.004
count = 0
chi_len = 15
Delta=8000.
kappa=1.
eps0=0.3
theta=100.
Omega=380. #DEFAULT


#print([np.min(U0_MD), np.max(U0_MD), np.min(Uf_MD), np.max(Uf_MD)])
#sys.exit()
"""
# Variation of Kappa
kappas = [0., 0.5, 1.0]

# Figure to store stress-strain evolution with kappa variation
f1, ax1 = plt.subplots(1,1)

for i, kappa in enumerate(kappas):
    if os.path.exists('sct_q.out'):
        command = ['rm', '-r', 'sct_q.out']
        subprocess.run(command)

    command = ['./shear_energy', 'qs', PE_0,
      '{}'.format(beta), 
      '{}'.format(u0),
      'chi_len', '{}'.format(chi_len),
      'kappa', '{}'.format(kappa),
      'c0', '{}'.format(c0),
      'eps0', '{}'.format(eps0),
      'Omega', '{}'.format(Omega),
      'theta', '{}'.format(theta),
      'Delta', '{}'.format(Delta),
      'chi_inf', '{}'.format(T_inf),
      'sy', '{}'.format(sy),
      'mu', '{}'.format(mu),
      'K', '{}'.format(K),
      'rho', '{}'.format(rho)]
    
    fig = plt.figure(figsize=(11,8.5))
    textStr = '\n'.join((
           r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
           r'$u_0$ = %.4f eV'%(u0, ),
           r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
           r'$\kappa$ = %.3f --' % (kappa, ),
           r'$c_0$ = %.3f --' % (c0, ),
           r'$\epsilon_0$ = %.3f --' % (eps0, ),
           r'$\Omega$ = %d $\AA^3$' % (Omega, ),
           r'$\theta$ = %d K' % (theta, ),
           r'$\Delta$ = %d K' % (Delta, ),
           r'$T_\infty$ = %d K' % (T_inf, ),
           r'$\sigma_y$ = %.2f GPa' % (sy, ),
           r'$\mu$ = %.2f GPa' % (mu, ),
           r'K = %.2f GPa' % (K, ),
           r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))

    ax = plt.subplot(334)
    ax.axis('off')
    ax.text(0, 0.5, textStr, verticalalignment='center')

    try:
      subprocess.run(command, timeout=300)
  
      if os.path.isfile('sct_q.out/Exy.100'):

        fig = plt.figure(figsize=(11,8.5))
        textStr = '\n'.join((
          r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
          r'$u_0$ = %.4f eV'%(u0, ),
          r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
          r'$\kappa$ = %.3f --' % (kappa, ),
          r'$c_0$ = %.3f --' % (c0, ),
          r'$\epsilon_0$ = %.3f --' % (eps0, ),
          r'$\Omega$ = %d $\AA^3$' % (Omega, ),
          r'$\theta$ = %d K' % (theta, ),
          r'$\Delta$ = %d K' % (Delta, ),
          r'$T_\infty$ = %d K' % (T_inf, ),
          r'$\sigma_y$ = %.2f GPa' % (sy, ),
          r'$\mu$ = %.2f GPa' % (mu, ),
          r'K = %.2f GPa' % (K, ),
          r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))



        ax = plt.subplot(334)
        ax.axis('off')
        ax.text(0, 0.5, textStr, verticalalignment='center')

        for t, it in enumerate(time):
          TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
          ny = int(TAU[0]) + 1
          nx = int(len(TAU)/ny)
          TAU = TAU.reshape(nx, ny)
          TAU = TAU[1:, 1:]
          
          tau[it] = np.mean(TAU)


        plt.subplot(332)
        plt.title(r'$\tau(\gamma)$')
        plt.plot(strain, tau, 'k-', label='CM')
        plt.plot(strain_md, tau_md, 'k:', label='MD')
        plt.xlabel(r'$\gamma$ [--]')
        plt.ylabel(r'$\tau$ [Gpa]')
        plt.legend()

        # Initial effective temperature
        T0 = np.fromfile(fileDir + 'tem.0', dtype=np.float32)
        ny = int(T0[0]) + 1
        nx = int(len(T0)/ny)
        T0 = T0.reshape(nx, ny)
        T0 = T0[1:, 1:]

        plt.subplot(333)
        plt.title(r'$T_{0}$')
        plt.contourf(T0, cmap='Greys')
        plt.colorbar()
        plt.axis('off')

        # Final CM effective temperature
        Tf_CM = np.fromfile(fileDir + 'tem.100', dtype=np.float32)
        Tf_CM = Tf_CM.reshape(nx, ny)
        Tf_CM = Tf_CM[1:, 1:]
        
        plt.subplot(335)
        plt.title(r'${T_{f}}^{CM}$')
        plt.contourf(Tf_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final CM plastic strain
        Exy_CM = np.fromfile(fileDir + 'Exy.100', dtype=np.float32)
        Exy_CM = Exy_CM.reshape(nx, ny)
        Exy_CM = 2.*Exy_CM[1:, 1:]

        plt.subplot(336)
        plt.title(r'${\Gamma_{f}}^{CM}$')
        plt.contourf(Exy_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD effective temperature
        Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
        Tf_MD = beta*(Uf_MD-u0)*TZ
        plt.subplot(338)
        plt.title(r'${T_{f}}^{MD}$')
        plt.contourf(Tf_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD plastic strain
        Exy_MD = np.loadtxt(refDir + 'Exy.MD.100', skiprows=1)
        plt.subplot(339)
        plt.title(r'${\Gamma_{f}}^{MD}$')
        plt.contourf(Exy_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Plot Stress-Strain
        if i == 0:
            ax1.plot(strain_md, tau_md, label='MD')

        ax1.plot(strain, tau, label=kappa)
        ax1.set_xlabel('Strain')
        ax1.set_ylabel('Stress')
        ax1.legend(title=r'$\kappa$')
        f1.savefig('kappa--stress_v_strain.png')        

        # Output data to screen
        textStr = '\n'.join((
            'c_0 = %.3f --' % (c0, ),
            'Omega = %d A**3' % (Omega, ),
            'Delta = %d K' % (Delta, ),
            'max(tau) = %.3f' % (np.amax(tau), ),
            'tau_ss = %.3f' % (tau[-1], ),
            'max(T_eff) = %.3f' % (np.amax(Tf_CM), ),
            'max(tau) = %.3f\n\n' % (np.amax(Exy_CM, ))))
        
        print(textStr)

        fig.tight_layout()

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1

    except subprocess.TimeoutExpired:

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1
      continue
"""

"""
# Variation of chi_len
chi_lens = [1., 10., 20.]
chi_lens = np.arange(1., 20., 2)
kappa - 1.
# Figure to store stress-strain evolution with kappa variation
f1, ax1 = plt.subplots(1,1)

for i, chi_len in enumerate(chi_lens):
    if os.path.exists('sct_q.out'):
        command = ['rm', '-r', 'sct_q.out']
        subprocess.run(command)

    command = ['./shear_energy', 'qs', PE_0,
      '{}'.format(beta), 
      '{}'.format(u0),
      'chi_len', '{}'.format(chi_len),
      'kappa', '{}'.format(kappa),
      'c0', '{}'.format(c0),
      'eps0', '{}'.format(eps0),
      'Omega', '{}'.format(Omega),
      'theta', '{}'.format(theta),
      'Delta', '{}'.format(Delta),
      'chi_inf', '{}'.format(T_inf),
      'sy', '{}'.format(sy),
      'mu', '{}'.format(mu),
      'K', '{}'.format(K),
      'rho', '{}'.format(rho)]
    
    fig = plt.figure(figsize=(11,8.5))
    textStr = '\n'.join((
           r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
           r'$u_0$ = %.4f eV'%(u0, ),
           r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
           r'$\kappa$ = %.3f --' % (kappa, ),
           r'$c_0$ = %.3f --' % (c0, ),
           r'$\epsilon_0$ = %.3f --' % (eps0, ),
           r'$\Omega$ = %d $\AA^3$' % (Omega, ),
           r'$\theta$ = %d K' % (theta, ),
           r'$\Delta$ = %d K' % (Delta, ),
           r'$T_\infty$ = %d K' % (T_inf, ),
           r'$\sigma_y$ = %.2f GPa' % (sy, ),
           r'$\mu$ = %.2f GPa' % (mu, ),
           r'K = %.2f GPa' % (K, ),
           r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))

    ax = plt.subplot(334)
    ax.axis('off')
    ax.text(0, 0.5, textStr, verticalalignment='center')

    try:
      subprocess.run(command, timeout=300)
  
      if os.path.isfile('sct_q.out/Exy.100'):

        fig = plt.figure(figsize=(11,8.5))
        textStr = '\n'.join((
          r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
          r'$u_0$ = %.4f eV'%(u0, ),
          r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
          r'$\kappa$ = %.3f --' % (kappa, ),
          r'$c_0$ = %.3f --' % (c0, ),
          r'$\epsilon_0$ = %.3f --' % (eps0, ),
          r'$\Omega$ = %d $\AA^3$' % (Omega, ),
          r'$\theta$ = %d K' % (theta, ),
          r'$\Delta$ = %d K' % (Delta, ),
          r'$T_\infty$ = %d K' % (T_inf, ),
          r'$\sigma_y$ = %.2f GPa' % (sy, ),
          r'$\mu$ = %.2f GPa' % (mu, ),
          r'K = %.2f GPa' % (K, ),
          r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))



        ax = plt.subplot(334)
        ax.axis('off')
        ax.text(0, 0.5, textStr, verticalalignment='center')

        for t, it in enumerate(time):
          TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
          ny = int(TAU[0]) + 1
          nx = int(len(TAU)/ny)
          TAU = TAU.reshape(nx, ny)
          TAU = TAU[1:, 1:]
          
          tau[it] = np.mean(TAU)


        plt.subplot(332)
        plt.title(r'$\tau(\gamma)$')
        plt.plot(strain, tau, 'k-', label='CM')
        plt.plot(strain_md, tau_md, 'k:', label='MD')
        plt.xlabel(r'$\gamma$ [--]')
        plt.ylabel(r'$\tau$ [Gpa]')
        plt.legend()

        # Initial effective temperature
        T0 = np.fromfile(fileDir + 'tem.0', dtype=np.float32)
        ny = int(T0[0]) + 1
        nx = int(len(T0)/ny)
        T0 = T0.reshape(nx, ny)
        T0 = T0[1:, 1:]

        plt.subplot(333)
        plt.title(r'$T_{0}$')
        plt.contourf(T0, cmap='Greys')
        plt.colorbar()
        plt.axis('off')

        # Final CM effective temperature
        Tf_CM = np.fromfile(fileDir + 'tem.100', dtype=np.float32)
        Tf_CM = Tf_CM.reshape(nx, ny)
        Tf_CM = Tf_CM[1:, 1:]
        
        plt.subplot(335)
        plt.title(r'${T_{f}}^{CM}$')
        plt.contourf(Tf_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final CM plastic strain
        Exy_CM = np.fromfile(fileDir + 'Exy.100', dtype=np.float32)
        Exy_CM = Exy_CM.reshape(nx, ny)
        Exy_CM = 2.*Exy_CM[1:, 1:]

        plt.subplot(336)
        plt.title(r'${\Gamma_{f}}^{CM}$')
        plt.contourf(Exy_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD effective temperature
        Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
        Tf_MD = beta*(Uf_MD-u0)*TZ
        plt.subplot(338)
        plt.title(r'${T_{f}}^{MD}$')
        plt.contourf(Tf_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD plastic strain
        Exy_MD = np.loadtxt(refDir + 'Exy.MD.100', skiprows=1)
        plt.subplot(339)
        plt.title(r'${\Gamma_{f}}^{MD}$')
        plt.contourf(Exy_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Plot Stress-Strain
        if i == 0:
            ax1.plot(strain_md, tau_md, label='MD')

        ax1.plot(strain, tau, label=chi_len)
        ax1.set_xlabel('Strain')
        ax1.set_ylabel('Stress')
        ax1.legend(title=r'$l_{\chi}$')
        f1.savefig('chi_len--stress_v_strain.png')        

        # Output data to screen
        textStr = '\n'.join((
            'c_0 = %.3f --' % (c0, ),
            'Omega = %d A**3' % (Omega, ),
            'Delta = %d K' % (Delta, ),
            'max(tau) = %.3f' % (np.amax(tau), ),
            'tau_ss = %.3f' % (tau[-1], ),
            'max(T_eff) = %.3f' % (np.amax(Tf_CM), ),
            'max(tau) = %.3f\n\n' % (np.amax(Exy_CM, ))))
        
        print(textStr)

        fig.tight_layout()

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1

    except subprocess.TimeoutExpired:

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1
      continue
"""

"""
# Variation of c0
c0s = np.arange(0.1, 0.6, 0.1)
# Figure to store stress-strain evolution with kappa variation
f1, ax1 = plt.subplots(1,1)

for i, c0 in enumerate(c0s):
    if os.path.exists('sct_q.out'):
        command = ['rm', '-r', 'sct_q.out']
        subprocess.run(command)

    command = ['./shear_energy', 'qs', PE_0,
      '{}'.format(beta), 
      '{}'.format(u0),
      'chi_len', '{}'.format(chi_len),
      'kappa', '{}'.format(kappa),
      'c0', '{}'.format(c0),
      'eps0', '{}'.format(eps0),
      'Omega', '{}'.format(Omega),
      'theta', '{}'.format(theta),
      'Delta', '{}'.format(Delta),
      'chi_inf', '{}'.format(T_inf),
      'sy', '{}'.format(sy),
      'mu', '{}'.format(mu),
      'K', '{}'.format(K),
      'rho', '{}'.format(rho)]
    
    fig = plt.figure(figsize=(11,8.5))
    textStr = '\n'.join((
           r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
           r'$u_0$ = %.4f eV'%(u0, ),
           r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
           r'$\kappa$ = %.3f --' % (kappa, ),
           r'$c_0$ = %.3f --' % (c0, ),
           r'$\epsilon_0$ = %.3f --' % (eps0, ),
           r'$\Omega$ = %d $\AA^3$' % (Omega, ),
           r'$\theta$ = %d K' % (theta, ),
           r'$\Delta$ = %d K' % (Delta, ),
           r'$T_\infty$ = %d K' % (T_inf, ),
           r'$\sigma_y$ = %.2f GPa' % (sy, ),
           r'$\mu$ = %.2f GPa' % (mu, ),
           r'K = %.2f GPa' % (K, ),
           r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))

    ax = plt.subplot(334)
    ax.axis('off')
    ax.text(0, 0.5, textStr, verticalalignment='center')

    try:
      subprocess.run(command, timeout=300)
  
      if os.path.isfile('sct_q.out/Exy.100'):

        fig = plt.figure(figsize=(11,8.5))
        textStr = '\n'.join((
          r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
          r'$u_0$ = %.4f eV'%(u0, ),
          r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
          r'$\kappa$ = %.3f --' % (kappa, ),
          r'$c_0$ = %.3f --' % (c0, ),
          r'$\epsilon_0$ = %.3f --' % (eps0, ),
          r'$\Omega$ = %d $\AA^3$' % (Omega, ),
          r'$\theta$ = %d K' % (theta, ),
          r'$\Delta$ = %d K' % (Delta, ),
          r'$T_\infty$ = %d K' % (T_inf, ),
          r'$\sigma_y$ = %.2f GPa' % (sy, ),
          r'$\mu$ = %.2f GPa' % (mu, ),
          r'K = %.2f GPa' % (K, ),
          r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))



        ax = plt.subplot(334)
        ax.axis('off')
        ax.text(0, 0.5, textStr, verticalalignment='center')

        for t, it in enumerate(time):
          TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
          ny = int(TAU[0]) + 1
          nx = int(len(TAU)/ny)
          TAU = TAU.reshape(nx, ny)
          TAU = TAU[1:, 1:]
          
          tau[it] = np.mean(TAU)


        plt.subplot(332)
        plt.title(r'$\tau(\gamma)$')
        plt.plot(strain, tau, 'k-', label='CM')
        plt.plot(strain_md, tau_md, 'k:', label='MD')
        plt.xlabel(r'$\gamma$ [--]')
        plt.ylabel(r'$\tau$ [Gpa]')
        plt.legend()

        # Initial effective temperature
        T0 = np.fromfile(fileDir + 'tem.0', dtype=np.float32)
        ny = int(T0[0]) + 1
        nx = int(len(T0)/ny)
        T0 = T0.reshape(nx, ny)
        T0 = T0[1:, 1:]

        plt.subplot(333)
        plt.title(r'$T_{0}$')
        plt.contourf(T0, cmap='Greys')
        plt.colorbar()
        plt.axis('off')

        # Final CM effective temperature
        Tf_CM = np.fromfile(fileDir + 'tem.100', dtype=np.float32)
        Tf_CM = Tf_CM.reshape(nx, ny)
        Tf_CM = Tf_CM[1:, 1:]
        
        plt.subplot(335)
        plt.title(r'${T_{f}}^{CM}$')
        plt.contourf(Tf_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final CM plastic strain
        Exy_CM = np.fromfile(fileDir + 'Exy.100', dtype=np.float32)
        Exy_CM = Exy_CM.reshape(nx, ny)
        Exy_CM = 2.*Exy_CM[1:, 1:]

        plt.subplot(336)
        plt.title(r'${\Gamma_{f}}^{CM}$')
        plt.contourf(Exy_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD effective temperature
        Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
        Tf_MD = beta*(Uf_MD-u0)*TZ
        plt.subplot(338)
        plt.title(r'${T_{f}}^{MD}$')
        plt.contourf(Tf_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD plastic strain
        Exy_MD = np.loadtxt(refDir + 'Exy.MD.100', skiprows=1)
        plt.subplot(339)
        plt.title(r'${\Gamma_{f}}^{MD}$')
        plt.contourf(Exy_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Plot Stress-Strain
        if i == 0:
            ax1.plot(strain_md, tau_md, label='MD')

        ax1.plot(strain, tau, label=c0)
        ax1.set_xlabel('Strain')
        ax1.set_ylabel('Stress')
        ax1.legend(title=r'$c_0$')
        f1.savefig('c0--stress_v_strain.png')        

        # Output data to screen
        textStr = '\n'.join((
            'c_0 = %.3f --' % (c0, ),
            'Omega = %d A**3' % (Omega, ),
            'Delta = %d K' % (Delta, ),
            'max(tau) = %.3f' % (np.amax(tau), ),
            'tau_ss = %.3f' % (tau[-1], ),
            'max(T_eff) = %.3f' % (np.amax(Tf_CM), ),
            'max(tau) = %.3f\n\n' % (np.amax(Exy_CM, ))))
        
        print(textStr)

        fig.tight_layout()

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1

    except subprocess.TimeoutExpired:

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1
      continue
"""

"""

# Variation of eps0
eps0s = np.arange(0.1, 0.6, 0.1)
# Figure to store stress-strain evolution with kappa variation
f1, ax1 = plt.subplots(1,1)

for i, eps0 in enumerate(eps0s):
    if os.path.exists('sct_q.out'):
        command = ['rm', '-r', 'sct_q.out']
        subprocess.run(command)

    command = ['./shear_energy', 'qs', PE_0,
      '{}'.format(beta), 
      '{}'.format(u0),
      'chi_len', '{}'.format(chi_len),
      'kappa', '{}'.format(kappa),
      'c0', '{}'.format(c0),
      'eps0', '{}'.format(eps0),
      'Omega', '{}'.format(Omega),
      'theta', '{}'.format(theta),
      'Delta', '{}'.format(Delta),
      'chi_inf', '{}'.format(T_inf),
      'sy', '{}'.format(sy),
      'mu', '{}'.format(mu),
      'K', '{}'.format(K),
      'rho', '{}'.format(rho)]
    
    fig = plt.figure(figsize=(11,8.5))
    textStr = '\n'.join((
           r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
           r'$u_0$ = %.4f eV'%(u0, ),
           r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
           r'$\kappa$ = %.3f --' % (kappa, ),
           r'$c_0$ = %.3f --' % (c0, ),
           r'$\epsilon_0$ = %.3f --' % (eps0, ),
           r'$\Omega$ = %d $\AA^3$' % (Omega, ),
           r'$\theta$ = %d K' % (theta, ),
           r'$\Delta$ = %d K' % (Delta, ),
           r'$T_\infty$ = %d K' % (T_inf, ),
           r'$\sigma_y$ = %.2f GPa' % (sy, ),
           r'$\mu$ = %.2f GPa' % (mu, ),
           r'K = %.2f GPa' % (K, ),
           r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))

    ax = plt.subplot(334)
    ax.axis('off')
    ax.text(0, 0.5, textStr, verticalalignment='center')

    try:
      subprocess.run(command, timeout=300)
  
      if os.path.isfile('sct_q.out/Exy.100'):

        fig = plt.figure(figsize=(11,8.5))
        textStr = '\n'.join((
          r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
          r'$u_0$ = %.4f eV'%(u0, ),
          r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
          r'$\kappa$ = %.3f --' % (kappa, ),
          r'$c_0$ = %.3f --' % (c0, ),
          r'$\epsilon_0$ = %.3f --' % (eps0, ),
          r'$\Omega$ = %d $\AA^3$' % (Omega, ),
          r'$\theta$ = %d K' % (theta, ),
          r'$\Delta$ = %d K' % (Delta, ),
          r'$T_\infty$ = %d K' % (T_inf, ),
          r'$\sigma_y$ = %.2f GPa' % (sy, ),
          r'$\mu$ = %.2f GPa' % (mu, ),
          r'K = %.2f GPa' % (K, ),
          r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))



        ax = plt.subplot(334)
        ax.axis('off')
        ax.text(0, 0.5, textStr, verticalalignment='center')

        for t, it in enumerate(time):
          TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
          ny = int(TAU[0]) + 1
          nx = int(len(TAU)/ny)
          TAU = TAU.reshape(nx, ny)
          TAU = TAU[1:, 1:]
          
          tau[it] = np.mean(TAU)


        plt.subplot(332)
        plt.title(r'$\tau(\gamma)$')
        plt.plot(strain, tau, 'k-', label='CM')
        plt.plot(strain_md, tau_md, 'k:', label='MD')
        plt.xlabel(r'$\gamma$ [--]')
        plt.ylabel(r'$\tau$ [Gpa]')
        plt.legend()

        # Initial effective temperature
        T0 = np.fromfile(fileDir + 'tem.0', dtype=np.float32)
        ny = int(T0[0]) + 1
        nx = int(len(T0)/ny)
        T0 = T0.reshape(nx, ny)
        T0 = T0[1:, 1:]

        plt.subplot(333)
        plt.title(r'$T_{0}$')
        plt.contourf(T0, cmap='Greys')
        plt.colorbar()
        plt.axis('off')

        # Final CM effective temperature
        Tf_CM = np.fromfile(fileDir + 'tem.100', dtype=np.float32)
        Tf_CM = Tf_CM.reshape(nx, ny)
        Tf_CM = Tf_CM[1:, 1:]
        
        plt.subplot(335)
        plt.title(r'${T_{f}}^{CM}$')
        plt.contourf(Tf_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final CM plastic strain
        Exy_CM = np.fromfile(fileDir + 'Exy.100', dtype=np.float32)
        Exy_CM = Exy_CM.reshape(nx, ny)
        Exy_CM = 2.*Exy_CM[1:, 1:]

        plt.subplot(336)
        plt.title(r'${\Gamma_{f}}^{CM}$')
        plt.contourf(Exy_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD effective temperature
        Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
        Tf_MD = beta*(Uf_MD-u0)*TZ
        plt.subplot(338)
        plt.title(r'${T_{f}}^{MD}$')
        plt.contourf(Tf_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD plastic strain
        Exy_MD = np.loadtxt(refDir + 'Exy.MD.100', skiprows=1)
        plt.subplot(339)
        plt.title(r'${\Gamma_{f}}^{MD}$')
        plt.contourf(Exy_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Plot Stress-Strain
        if i == 0:
            ax1.plot(strain_md, tau_md, label='MD')

        ax1.plot(strain, tau, label=eps0)
        ax1.set_xlabel('Strain')
        ax1.set_ylabel('Stress')
        ax1.legend(title=r'$\epsilon_0$')
        f1.savefig('eps0--stress_v_strain.png')        

        # Output data to screen
        textStr = '\n'.join((
            'c_0 = %.3f --' % (c0, ),
            'Omega = %d A**3' % (Omega, ),
            'Delta = %d K' % (Delta, ),
            'max(tau) = %.3f' % (np.amax(tau), ),
            'tau_ss = %.3f' % (tau[-1], ),
            'max(T_eff) = %.3f' % (np.amax(Tf_CM), ),
            'max(tau) = %.3f\n\n' % (np.amax(Exy_CM, ))))
        
        print(textStr)

        fig.tight_layout()

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1

    except subprocess.TimeoutExpired:

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1
      continue
"""

"""
# Variation of Omega
Omegas = np.arange(360, 500, 20)
# Figure to store stress-strain evolution with kappa variation
f1, ax1 = plt.subplots(1,1)

for i, Omega in enumerate(Omegas):
    if os.path.exists('sct_q.out'):
        command = ['rm', '-r', 'sct_q.out']
        subprocess.run(command)

    command = ['./shear_energy', 'qs', PE_0,
      '{}'.format(beta), 
      '{}'.format(u0),
      'chi_len', '{}'.format(chi_len),
      'kappa', '{}'.format(kappa),
      'c0', '{}'.format(c0),
      'eps0', '{}'.format(eps0),
      'Omega', '{}'.format(Omega),
      'theta', '{}'.format(theta),
      'Delta', '{}'.format(Delta),
      'chi_inf', '{}'.format(T_inf),
      'sy', '{}'.format(sy),
      'mu', '{}'.format(mu),
      'K', '{}'.format(K),
      'rho', '{}'.format(rho)]
    
    fig = plt.figure(figsize=(11,8.5))
    textStr = '\n'.join((
           r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
           r'$u_0$ = %.4f eV'%(u0, ),
           r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
           r'$\kappa$ = %.3f --' % (kappa, ),
           r'$c_0$ = %.3f --' % (c0, ),
           r'$\epsilon_0$ = %.3f --' % (eps0, ),
           r'$\Omega$ = %d $\AA^3$' % (Omega, ),
           r'$\theta$ = %d K' % (theta, ),
           r'$\Delta$ = %d K' % (Delta, ),
           r'$T_\infty$ = %d K' % (T_inf, ),
           r'$\sigma_y$ = %.2f GPa' % (sy, ),
           r'$\mu$ = %.2f GPa' % (mu, ),
           r'K = %.2f GPa' % (K, ),
           r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))

    ax = plt.subplot(334)
    ax.axis('off')
    ax.text(0, 0.5, textStr, verticalalignment='center')

    try:
      subprocess.run(command, timeout=300)
  
      if os.path.isfile('sct_q.out/Exy.100'):

        fig = plt.figure(figsize=(11,8.5))
        textStr = '\n'.join((
          r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
          r'$u_0$ = %.4f eV'%(u0, ),
          r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
          r'$\kappa$ = %.3f --' % (kappa, ),
          r'$c_0$ = %.3f --' % (c0, ),
          r'$\epsilon_0$ = %.3f --' % (eps0, ),
          r'$\Omega$ = %d $\AA^3$' % (Omega, ),
          r'$\theta$ = %d K' % (theta, ),
          r'$\Delta$ = %d K' % (Delta, ),
          r'$T_\infty$ = %d K' % (T_inf, ),
          r'$\sigma_y$ = %.2f GPa' % (sy, ),
          r'$\mu$ = %.2f GPa' % (mu, ),
          r'K = %.2f GPa' % (K, ),
          r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))



        ax = plt.subplot(334)
        ax.axis('off')
        ax.text(0, 0.5, textStr, verticalalignment='center')

        for t, it in enumerate(time):
          TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
          ny = int(TAU[0]) + 1
          nx = int(len(TAU)/ny)
          TAU = TAU.reshape(nx, ny)
          TAU = TAU[1:, 1:]
          
          tau[it] = np.mean(TAU)


        plt.subplot(332)
        plt.title(r'$\tau(\gamma)$')
        plt.plot(strain, tau, 'k-', label='CM')
        plt.plot(strain_md, tau_md, 'k:', label='MD')
        plt.xlabel(r'$\gamma$ [--]')
        plt.ylabel(r'$\tau$ [Gpa]')
        plt.legend()

        # Initial effective temperature
        T0 = np.fromfile(fileDir + 'tem.0', dtype=np.float32)
        ny = int(T0[0]) + 1
        nx = int(len(T0)/ny)
        T0 = T0.reshape(nx, ny)
        T0 = T0[1:, 1:]

        plt.subplot(333)
        plt.title(r'$T_{0}$')
        plt.contourf(T0, cmap='Greys')
        plt.colorbar()
        plt.axis('off')

        # Final CM effective temperature
        Tf_CM = np.fromfile(fileDir + 'tem.100', dtype=np.float32)
        Tf_CM = Tf_CM.reshape(nx, ny)
        Tf_CM = Tf_CM[1:, 1:]
        
        plt.subplot(335)
        plt.title(r'${T_{f}}^{CM}$')
        plt.contourf(Tf_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final CM plastic strain
        Exy_CM = np.fromfile(fileDir + 'Exy.100', dtype=np.float32)
        Exy_CM = Exy_CM.reshape(nx, ny)
        Exy_CM = 2.*Exy_CM[1:, 1:]

        plt.subplot(336)
        plt.title(r'${\Gamma_{f}}^{CM}$')
        plt.contourf(Exy_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD effective temperature
        Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
        Tf_MD = beta*(Uf_MD-u0)*TZ
        plt.subplot(338)
        plt.title(r'${T_{f}}^{MD}$')
        plt.contourf(Tf_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD plastic strain
        Exy_MD = np.loadtxt(refDir + 'Exy.MD.100', skiprows=1)
        plt.subplot(339)
        plt.title(r'${\Gamma_{f}}^{MD}$')
        plt.contourf(Exy_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Plot Stress-Strain
        if i == 0:
            ax1.plot(strain_md, tau_md, label='MD')

        ax1.plot(strain, tau, label=Omega)
        ax1.set_xlabel('Strain')
        ax1.set_ylabel('Stress')
        ax1.legend(title=r'$\Omega$')
        f1.savefig('Omega--stress_v_strain.png')        

        # Output data to screen
        textStr = '\n'.join((
            'c_0 = %.3f --' % (c0, ),
            'Omega = %d A**3' % (Omega, ),
            'Delta = %d K' % (Delta, ),
            'max(tau) = %.3f' % (np.amax(tau), ),
            'tau_ss = %.3f' % (tau[-1], ),
            'max(T_eff) = %.3f' % (np.amax(Tf_CM), ),
            'max(tau) = %.3f\n\n' % (np.amax(Exy_CM, ))))
        
        print(textStr)

        fig.tight_layout()

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1

    except subprocess.TimeoutExpired:

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1
      continue
"""
"""
# Variation of Delta
Deltas = np.arange(7800, 9000, 200)
# Figure to store stress-strain evolution with kappa variation
f1, ax1 = plt.subplots(1,1)

for i, Delta in enumerate(Deltas):
    if os.path.exists('sct_q.out'):
        command = ['rm', '-r', 'sct_q.out']
        subprocess.run(command)

    command = ['./shear_energy', 'qs', PE_0,
      '{}'.format(beta), 
      '{}'.format(u0),
      'chi_len', '{}'.format(chi_len),
      'kappa', '{}'.format(kappa),
      'c0', '{}'.format(c0),
      'eps0', '{}'.format(eps0),
      'Omega', '{}'.format(Omega),
      'theta', '{}'.format(theta),
      'Delta', '{}'.format(Delta),
      'chi_inf', '{}'.format(T_inf),
      'sy', '{}'.format(sy),
      'mu', '{}'.format(mu),
      'K', '{}'.format(K),
      'rho', '{}'.format(rho)]
    
    fig = plt.figure(figsize=(11,8.5))
    textStr = '\n'.join((
           r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
           r'$u_0$ = %.4f eV'%(u0, ),
           r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
           r'$\kappa$ = %.3f --' % (kappa, ),
           r'$c_0$ = %.3f --' % (c0, ),
           r'$\epsilon_0$ = %.3f --' % (eps0, ),
           r'$\Omega$ = %d $\AA^3$' % (Omega, ),
           r'$\theta$ = %d K' % (theta, ),
           r'$\Delta$ = %d K' % (Delta, ),
           r'$T_\infty$ = %d K' % (T_inf, ),
           r'$\sigma_y$ = %.2f GPa' % (sy, ),
           r'$\mu$ = %.2f GPa' % (mu, ),
           r'K = %.2f GPa' % (K, ),
           r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))

    ax = plt.subplot(334)
    ax.axis('off')
    ax.text(0, 0.5, textStr, verticalalignment='center')

    try:
      subprocess.run(command, timeout=300)
  
      if os.path.isfile('sct_q.out/Exy.100'):

        fig = plt.figure(figsize=(11,8.5))
        textStr = '\n'.join((
          r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
          r'$u_0$ = %.4f eV'%(u0, ),
          r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
          r'$\kappa$ = %.3f --' % (kappa, ),
          r'$c_0$ = %.3f --' % (c0, ),
          r'$\epsilon_0$ = %.3f --' % (eps0, ),
          r'$\Omega$ = %d $\AA^3$' % (Omega, ),
          r'$\theta$ = %d K' % (theta, ),
          r'$\Delta$ = %d K' % (Delta, ),
          r'$T_\infty$ = %d K' % (T_inf, ),
          r'$\sigma_y$ = %.2f GPa' % (sy, ),
          r'$\mu$ = %.2f GPa' % (mu, ),
          r'K = %.2f GPa' % (K, ),
          r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))



        ax = plt.subplot(334)
        ax.axis('off')
        ax.text(0, 0.5, textStr, verticalalignment='center')

        for t, it in enumerate(time):
          TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
          ny = int(TAU[0]) + 1
          nx = int(len(TAU)/ny)
          TAU = TAU.reshape(nx, ny)
          TAU = TAU[1:, 1:]
          
          tau[it] = np.mean(TAU)


        plt.subplot(332)
        plt.title(r'$\tau(\gamma)$')
        plt.plot(strain, tau, 'k-', label='CM')
        plt.plot(strain_md, tau_md, 'k:', label='MD')
        plt.xlabel(r'$\gamma$ [--]')
        plt.ylabel(r'$\tau$ [Gpa]')
        plt.legend()

        # Initial effective temperature
        T0 = np.fromfile(fileDir + 'tem.0', dtype=np.float32)
        ny = int(T0[0]) + 1
        nx = int(len(T0)/ny)
        T0 = T0.reshape(nx, ny)
        T0 = T0[1:, 1:]

        plt.subplot(333)
        plt.title(r'$T_{0}$')
        plt.contourf(T0, cmap='Greys')
        plt.colorbar()
        plt.axis('off')

        # Final CM effective temperature
        Tf_CM = np.fromfile(fileDir + 'tem.100', dtype=np.float32)
        Tf_CM = Tf_CM.reshape(nx, ny)
        Tf_CM = Tf_CM[1:, 1:]
        
        plt.subplot(335)
        plt.title(r'${T_{f}}^{CM}$')
        plt.contourf(Tf_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final CM plastic strain
        Exy_CM = np.fromfile(fileDir + 'Exy.100', dtype=np.float32)
        Exy_CM = Exy_CM.reshape(nx, ny)
        Exy_CM = 2.*Exy_CM[1:, 1:]

        plt.subplot(336)
        plt.title(r'${\Gamma_{f}}^{CM}$')
        plt.contourf(Exy_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD effective temperature
        Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
        Tf_MD = beta*(Uf_MD-u0)*TZ
        plt.subplot(338)
        plt.title(r'${T_{f}}^{MD}$')
        plt.contourf(Tf_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD plastic strain
        Exy_MD = np.loadtxt(refDir + 'Exy.MD.100', skiprows=1)
        plt.subplot(339)
        plt.title(r'${\Gamma_{f}}^{MD}$')
        plt.contourf(Exy_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Plot Stress-Strain
        if i == 0:
            ax1.plot(strain_md, tau_md, label='MD')

        ax1.plot(strain, tau, label=Delta)
        ax1.set_xlabel('Strain')
        ax1.set_ylabel('Stress')
        ax1.legend(title=r'$\Delta$')
        f1.savefig('Delta--stress_v_strain.png')        

        # Output data to screen
        textStr = '\n'.join((
            'c_0 = %.3f --' % (c0, ),
            'Omega = %d A**3' % (Omega, ),
            'Delta = %d K' % (Delta, ),
            'max(tau) = %.3f' % (np.amax(tau), ),
            'tau_ss = %.3f' % (tau[-1], ),
            'max(T_eff) = %.3f' % (np.amax(Tf_CM), ),
            'max(tau) = %.3f\n\n' % (np.amax(Exy_CM, ))))
        
        print(textStr)

        fig.tight_layout()

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1

    except subprocess.TimeoutExpired:

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1
      continue
"""

"""
# Variation of theta
thetas = np.arange(50, 120, 10)
# Figure to store stress-strain evolution with kappa variation
f1, ax1 = plt.subplots(1,1)

for i, theta in enumerate(thetas):
    if os.path.exists('sct_q.out'):
        command = ['rm', '-r', 'sct_q.out']
        subprocess.run(command)

    command = ['./shear_energy', 'qs', PE_0,
      '{}'.format(beta), 
      '{}'.format(u0),
      'chi_len', '{}'.format(chi_len),
      'kappa', '{}'.format(kappa),
      'c0', '{}'.format(c0),
      'eps0', '{}'.format(eps0),
      'Omega', '{}'.format(Omega),
      'theta', '{}'.format(theta),
      'Delta', '{}'.format(Delta),
      'chi_inf', '{}'.format(T_inf),
      'sy', '{}'.format(sy),
      'mu', '{}'.format(mu),
      'K', '{}'.format(K),
      'rho', '{}'.format(rho)]
    
    fig = plt.figure(figsize=(11,8.5))
    textStr = '\n'.join((
           r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
           r'$u_0$ = %.4f eV'%(u0, ),
           r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
           r'$\kappa$ = %.3f --' % (kappa, ),
           r'$c_0$ = %.3f --' % (c0, ),
           r'$\epsilon_0$ = %.3f --' % (eps0, ),
           r'$\Omega$ = %d $\AA^3$' % (Omega, ),
           r'$\theta$ = %d K' % (theta, ),
           r'$\Delta$ = %d K' % (Delta, ),
           r'$T_\infty$ = %d K' % (T_inf, ),
           r'$\sigma_y$ = %.2f GPa' % (sy, ),
           r'$\mu$ = %.2f GPa' % (mu, ),
           r'K = %.2f GPa' % (K, ),
           r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))

    ax = plt.subplot(334)
    ax.axis('off')
    ax.text(0, 0.5, textStr, verticalalignment='center')

    try:
      subprocess.run(command, timeout=300)
  
      if os.path.isfile('sct_q.out/Exy.100'):

        fig = plt.figure(figsize=(11,8.5))
        textStr = '\n'.join((
          r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
          r'$u_0$ = %.4f eV'%(u0, ),
          r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
          r'$\kappa$ = %.3f --' % (kappa, ),
          r'$c_0$ = %.3f --' % (c0, ),
          r'$\epsilon_0$ = %.3f --' % (eps0, ),
          r'$\Omega$ = %d $\AA^3$' % (Omega, ),
          r'$\theta$ = %d K' % (theta, ),
          r'$\Delta$ = %d K' % (Delta, ),
          r'$T_\infty$ = %d K' % (T_inf, ),
          r'$\sigma_y$ = %.2f GPa' % (sy, ),
          r'$\mu$ = %.2f GPa' % (mu, ),
          r'K = %.2f GPa' % (K, ),
          r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, )))



        ax = plt.subplot(334)
        ax.axis('off')
        ax.text(0, 0.5, textStr, verticalalignment='center')

        for t, it in enumerate(time):
          TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
          ny = int(TAU[0]) + 1
          nx = int(len(TAU)/ny)
          TAU = TAU.reshape(nx, ny)
          TAU = TAU[1:, 1:]
          
          tau[it] = np.mean(TAU)


        plt.subplot(332)
        plt.title(r'$\tau(\gamma)$')
        plt.plot(strain, tau, 'k-', label='CM')
        plt.plot(strain_md, tau_md, 'k:', label='MD')
        plt.xlabel(r'$\gamma$ [--]')
        plt.ylabel(r'$\tau$ [Gpa]')
        plt.legend()

        # Initial effective temperature
        T0 = np.fromfile(fileDir + 'tem.0', dtype=np.float32)
        ny = int(T0[0]) + 1
        nx = int(len(T0)/ny)
        T0 = T0.reshape(nx, ny)
        T0 = T0[1:, 1:]

        plt.subplot(333)
        plt.title(r'$T_{0}$')
        plt.contourf(T0, cmap='Greys')
        plt.colorbar()
        plt.axis('off')

        # Final CM effective temperature
        Tf_CM = np.fromfile(fileDir + 'tem.100', dtype=np.float32)
        Tf_CM = Tf_CM.reshape(nx, ny)
        Tf_CM = Tf_CM[1:, 1:]
        
        plt.subplot(335)
        plt.title(r'${T_{f}}^{CM}$')
        plt.contourf(Tf_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final CM plastic strain
        Exy_CM = np.fromfile(fileDir + 'Exy.100', dtype=np.float32)
        Exy_CM = Exy_CM.reshape(nx, ny)
        Exy_CM = 2.*Exy_CM[1:, 1:]

        plt.subplot(336)
        plt.title(r'${\Gamma_{f}}^{CM}$')
        plt.contourf(Exy_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD effective temperature
        Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
        Tf_MD = beta*(Uf_MD-u0)*TZ
        plt.subplot(338)
        plt.title(r'${T_{f}}^{MD}$')
        plt.contourf(Tf_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD plastic strain
        Exy_MD = np.loadtxt(refDir + 'Exy.MD.100', skiprows=1)
        plt.subplot(339)
        plt.title(r'${\Gamma_{f}}^{MD}$')
        plt.contourf(Exy_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Plot Stress-Strain
        if i == 0:
            ax1.plot(strain_md, tau_md, label='MD')

        ax1.plot(strain, tau, label=theta)
        ax1.set_xlabel('Strain')
        ax1.set_ylabel('Stress')
        ax1.legend(title=r'$\theta$')
        f1.savefig('theta--stress_v_strain.png')        

        # Output data to screen
        textStr = '\n'.join((
            'c_0 = %.3f --' % (c0, ),
            'Omega = %d A**3' % (Omega, ),
            'Delta = %d K' % (Delta, ),
            'max(tau) = %.3f' % (np.amax(tau), ),
            'tau_ss = %.3f' % (tau[-1], ),
            'max(T_eff) = %.3f' % (np.amax(Tf_CM), ),
            'max(tau) = %.3f\n\n' % (np.amax(Exy_CM, ))))
        
        print(textStr)

        fig.tight_layout()

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1

    except subprocess.TimeoutExpired:

      figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--eps0_{}--theta_{}--Omega_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, eps0, theta, Omega, Delta, sy, T_inf, rho, mu)
      plt.savefig(figTitle)
      plt.close()

      count += 1
      continue
"""

# Variation of u0
beta = 2.0
u0_min = np.amin([Uf_MD, U0_MD])
#print(u0_min)
#sys.exit()
u0s = np.arange(np.amin(Uf_MD), np.amin(Uf_MD)+0.005, 0.00001)
# Good for beta = 9.0
beta = 9.0
u0s = np.arange(-3.363, -3.361, 0.0002)

# Good for beta = 2.0
beta = 2.0
u0s = np.arange(-3.38, -3.36, 0.005)

# Good for beta = 9.5
beta = 9.5
u0s = np.arange(-3.362, -3.36, 0.0001)

# Good for beta = 5
beta = 5.0
u0s = np.arange(-3.38, -3.36, 0.005)
#u0s = np.arange(-3.363, -3.361, 0.00002)
#u0s = np.arange(-3.369, -3.361, 0.00002)

# Good for beta = 1
beta = 1.
u0s = np.arange(-3.39, -3.37, 0.0001)

# Good for beta = 9.5
beta = 9.5
u0s = np.arange(-3.39, -3.36, 0.0001)
## Good for beta = 10
#beta = 20.
#u0s = np.arange(-3.369, -3.361, 0.0001)

Tf_MD = beta*(Uf_MD-u0)*TZ
T_inf = np.amax(Tf_MD)

# Reference
#betas = np.arange(2,10,1).tolist()
beta=8
u0 = -3.361
#u0s = np.arange(-3.369,-3.360,0.001).tolist()
Deltas= np.arange(8000,9000,1000).tolist()
Omega=370
sy=0.85
theta=100
eps0=0.1
c0=1.5
rho=7234
K=143.23
mu=20
chi_len=1.937
Omegaeps0=111

# Figure to store stress-strain evolution with kappa variation
f1, ax1 = plt.subplots(1,1)

for i, Delta in enumerate(Deltas):

    Tf_MD = beta*(Uf_MD-u0)*TZ
    T_inf = np.amax(Tf_MD)

    if os.path.exists('sct_q.out'):
        command = ['rm', '-r', 'sct_q.out']
        subprocess.run(command)

    command = ['./shear_data', 'qs', PE_0,
      '{}'.format(beta), 
      '{}'.format(u0),
      'chi_len', '{}'.format(chi_len),
      'kappa', '{}'.format(kappa),
      'c0', '{}'.format(c0),
      'eps0', '{}'.format(eps0),
      'Omega', '{}'.format(Omega),
      'theta', '{}'.format(theta),
      'Delta', '{}'.format(Delta),
      'chi_inf', '{}'.format(T_inf),
      'sy', '{}'.format(sy),
      'mu', '{}'.format(mu),
      'K', '{}'.format(K),
      'rho', '{}'.format(rho),
      'Omegaeps0', '{}'.format(Omegaeps0)]
    
    fig = plt.figure(figsize=(11,8.5))
    textStr = '\n'.join((
           r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
           r'$u_0$ = %.4f eV'%(u0, ),
           r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
           r'$\kappa$ = %.3f --' % (kappa, ),
           r'$c_0$ = %.3f --' % (c0, ),
           r'$\epsilon_0$ = %.3f --' % (eps0, ),
           r'$\Omega$ = %d $\AA^3$' % (Omega, ),
           r'$\theta$ = %d K' % (theta, ),
           r'$\Delta$ = %d K' % (Delta, ),
           r'$T_\infty$ = %d K' % (T_inf, ),
           r'$\sigma_y$ = %.2f GPa' % (sy, ),
           r'$\mu$ = %.2f GPa' % (mu, ),
           r'K = %.2f GPa' % (K, ),
           r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, ),
           r'$Omegaeps0$ = %d $\AA^3$' % (Omegaeps0, )))

    ax = plt.subplot(334)
    ax.axis('off')
    ax.text(0, 0.5, textStr, verticalalignment='center')

    try:
      subprocess.run(command, timeout=30)
  
      if os.path.isfile('sct_q.out/Exy.100'):

        fig = plt.figure(figsize=(11,8.5))
        textStr = '\n'.join((
          r'$\beta$ = %.4f eV$^{-1}$'%(beta, ),
          r'$u_0$ = %.4f eV'%(u0, ),
          r'$l_\chi$ = %.3f $\AA$' % (chi_len, ), 
          r'$\kappa$ = %.3f --' % (kappa, ),
          r'$c_0$ = %.3f --' % (c0, ),
          r'$\epsilon_0$ = %.3f --' % (eps0, ),
          r'$\Omega$ = %d $\AA^3$' % (Omega, ),
          r'$\theta$ = %d K' % (theta, ),
          r'$\Delta$ = %d K' % (Delta, ),
          r'$T_\infty$ = %d K' % (T_inf, ),
          r'$\sigma_y$ = %.2f GPa' % (sy, ),
          r'$\mu$ = %.2f GPa' % (mu, ),
          r'K = %.2f GPa' % (K, ),
          r'$\rho$ = %d kg $\cdot m^{-3}$' % (rho, ),
          r'$Omegaeps0$ = %d $\AA^3$' % (Omegaeps0, )))



        ax = plt.subplot(334)
        ax.axis('off')
        ax.text(0, 0.5, textStr, verticalalignment='center')

        for t, it in enumerate(time):
          TAU = np.fromfile(fileDir + 'tau.' + str(t), dtype=np.float32)
          ny = int(TAU[0]) + 1
          nx = int(len(TAU)/ny)
          TAU = TAU.reshape(nx, ny)
          TAU = TAU[1:, 1:]
          
          tau[it] = np.mean(TAU)

          Q = np.fromfile(fileDir + 'q.' + str(t), dtype=np.float32)
          ny = int(Q[0]) + 1
          nx = int(len(Q)/ny)
          Q = Q.reshape(nx, ny)
          Q = Q[1:, 1:]
          
          q[it] = np.mean(Q)

          S = np.fromfile(fileDir + 's.' + str(t), dtype=np.float32)
          ny = int(S[0]) + 1
          nx = int(len(S)/ny)
          S = S.reshape(nx, ny)
          S = S[1:, 1:]
          
          s[it] = np.mean(S)

          dev[it] = np.sqrt(3*q[it]**2+s[it]**2+tau[it]**2)


        plt.subplot(332)
        plt.title(r'$\tau(\gamma)$')
        plt.plot(strain, tau, 'k-', label='CM')
        plt.plot(strain_md, tau_md, 'k:', label='MD')
        plt.xlabel(r'$\gamma$ [--]')
        plt.ylabel(r'$\tau$ [Gpa]')
        plt.legend()

        # Initial effective temperature
        T0 = np.fromfile(fileDir + 'tem.0', dtype=np.float32)
        ny = int(T0[0]) + 1
        nx = int(len(T0)/ny)
        T0 = T0.reshape(nx, ny)
        T0 = T0[1:, 1:]

        plt.subplot(333)
        plt.title(r'$T_{0}$')
        plt.contourf(T0, cmap='Greys')
        plt.colorbar()
        plt.axis('off')

        # Final CM effective temperature
        Tf_CM = np.fromfile(fileDir + 'tem.100', dtype=np.float32)
        Tf_CM = Tf_CM.reshape(nx, ny)
        Tf_CM = Tf_CM[1:, 1:]
        
        plt.subplot(335)
        plt.title(r'${T_{f}}^{CM}$')
        plt.contourf(Tf_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final CM plastic strain
        Exy_CM = np.fromfile(fileDir + 'Exy.100', dtype=np.float32)
        Exy_CM = Exy_CM.reshape(nx, ny)
        Exy_CM = 2.*Exy_CM[1:, 1:]

        plt.subplot(336)
        plt.title(r'${\Gamma_{f}}^{CM}$')
        plt.contourf(Exy_CM, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD effective temperature
        Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
        Tf_MD = beta*(Uf_MD-u0)*TZ
        plt.subplot(338)
        plt.title(r'${T_{f}}^{MD}$')
        plt.contourf(Tf_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Final MD plastic strain
        Exy_MD = np.loadtxt(refDir + 'Exy.MD.100', skiprows=1)
        plt.subplot(339)
        plt.title(r'${\Gamma_{f}}^{MD}$')
        plt.contourf(Exy_MD, cmap='Greys')
        plt.colorbar()
        plt.axis('off')
        
        # Plot Stress-Strain
        ax1.plot(strain_md, tau_md, label='MD')
        ax1.plot(strain, tau, label=u0)
        ax1.set_xlabel('Strain')
        ax1.set_ylabel('Stress')
        ax1.legend(title=r'$u0$')
        f1.savefig('u0--stress_v_strain.png')        

        # Output data to screen
        textStr = '\n'.join((
            'c_0 = %.3f --' % (c0, ),
            'Omega = %d A**3' % (Omega, ),
            'Delta = %d K' % (Delta, ),
            'max(tau) = %.3f' % (np.amax(tau), ),
            'tau_ss = %.3f' % (tau[-1], ),
            'max(T_eff) = %.3f' % (np.amax(Tf_CM), ),
            'max(tau) = %.3f\n\n' % (np.amax(Exy_CM, ))))
        
        print(textStr)

        fig.tight_layout()

        figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--Omegaeps0_{}--theta_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, Omegaeps0, theta, Delta, sy, T_inf, rho, mu)
        plt.savefig(figTitle)
        plt.close()

        count += 1

    except subprocess.TimeoutExpired:

        figTitle='Fields--beta_{:.4f}--u0_{:.4f}--chi_{}--kappa_{}--c0_{}--Omegaeps0_{}--theta_{}--Delta_{}--sy_{:.2f}--Tinf_{:.2f}--rho_{}--mu_{}.png'.format(beta, u0, chi_len, kappa, c0, Omegaeps0, theta, Delta, sy, T_inf, rho, mu)
        plt.savefig(figTitle)
        plt.close()

        count += 1

