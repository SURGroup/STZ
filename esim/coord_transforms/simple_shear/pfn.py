import numpy as np
import matplotlib.pyplot as plt
#plt.switch_backend('Qt5Agg')
#plt.ion()
import os, sys, fnmatch, subprocess

def STZmodel(samples):

    #plt.switch_backend('Agg')
    # Parameters and constants
    Ez = 1.809642
    TZ = 21000
    kB = 8.617330350e-5

    #BETA : (2 - 18eV)
    #U0 :   (np.amin(pe.MD.0) - 0.01)  —  (np.amax(pe.MD.0) - 0.001)
    #CHI_LEN :10 -100 Angstrom
    #KAPPA :10 -100 Angstrom
    #C0 :0.1 - 0.9
    #EPS0 :0.1 - 0.9
    #OMEGA :60 - 300 Cubic Angstroms
    #THETA :100K
    #DELTA :500 - 10000K
    #CHI_INF : BETA * (np.amax(pe.MD.100) - U0 ) * TZ, TZ = 21000
    #SY :0.85Gpa
    #MU :20GPa
    #K :143.23GPa
    #RHO :7234

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

    #refDir="/Users/dimitris/Google Drive/STZ/CuZr_Ref/"
    refDir="/Users/dimitrisgiovanis/Google Drive/STZ/CuZr_Ref/"
    fileDir = 'sct_q.out/'
    PE_0 = refDir+'pe.MD.0'

    time = np.arange(0, 101)
    tau = np.zeros(time.shape)
    s = np.zeros(time.shape)
    q = np.zeros(time.shape)
    tau_cm = np.zeros(time.shape)
    dev_cm = np.zeros(time.shape)

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
    #u0=u_max-0.0035
    #u0=u_max-0.004
    #print(u0)
    #count = 0

    beta = 4.7958
    u0 = -3.3625
    #chi_len =15.05
    chi_len = 10 #use 2 (increase this reduce the resolution)
    kappa = 0.06
    c0 = 0.414
    eps0 = 0.333
    Omega =349
    theta  = 97
    Delta = 5000#7948
    sy = 0.7
    mu = 20
    K = 143.23
    rho = 7234

    #rf_mean = 280 #KETSON
    #rf_std = 20 #KETSON
    print(samples)
    rf_mean = samples[0][0]
    rf_std = samples[0][1]
    
    Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
    Tf_MD = beta*(Uf_MD-u0)*TZ
    T_inf = np.amax(Tf_MD)

    if os.path.exists('sct_q.out'):
        command = ['rm', '-r', 'sct_q.out']
        subprocess.run(command)

    command = ['./shear_energy', 'qs', PE_0,
      '{}'.format(beta), 
      '{}'.format(u0),
      'chi_len', '{}'.format(chi_len),
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
      'rf_mean', '{}'.format(rf_mean),
      'rf_std', '{}'.format(rf_std)] #KETSON

    try:
      subprocess.run(command, timeout=300)

      if os.path.isfile('sct_q.out/Exy.100'):

        for t, it in enumerate(time):
          dev_ = np.fromfile(fileDir + 'dev.' + str(t), dtype=np.float32)
          ny = int(dev_[0]) + 1
          nx = int(len(dev_)/ny)
          dev_ = dev_.reshape(nx, ny)
          dev_ = dev_[1:, 1:]
          dev_cm[it] = np.mean(dev_)

        # Initial effective temperature
        T0 = np.fromfile(fileDir + 'tem.0', dtype=np.float32)
        ny = int(T0[0]) + 1
        nx = int(len(T0)/ny)
        T0 = T0.reshape(nx, ny)
        T0 = T0[1:, 1:]

        # Final CM effective temperature
        Tf_CM = np.fromfile(fileDir + 'tem.100', dtype=np.float32)
        Tf_CM = Tf_CM.reshape(nx, ny)
        Tf_CM = Tf_CM[1:, 1:]

        # Final CM plastic strain
        Exy_CM = np.fromfile(fileDir + 'Exy.100', dtype=np.float32)
        Exy_CM = Exy_CM.reshape(nx, ny)
        Exy_CM = 2.*Exy_CM[1:, 1:]

        # Final MD effective temperature
        Uf_MD = np.loadtxt(refDir + 'pe.MD.100', skiprows=1)
        Tf_MD = beta*(Uf_MD-u0)*TZ

        # Final MD plastic strain
        Exy_MD = np.loadtxt(refDir + 'Exy.MD.100', skiprows=1)

    except subprocess.TimeoutExpired:

           print('Timeout Expired')

    return Exy_CM
