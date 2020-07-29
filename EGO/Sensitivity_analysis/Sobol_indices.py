#!/usr/bin/env python

import numpy as np
import random
import shutil, sys, os, subprocess
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from scipy.stats import norm
import chaospy as cp
from UQpy.SampleMethods import LHS
from UQpy.Distributions import Uniform

sourceDir='/Users/katianakontolati/Projects/continuum_optimization/STZ/CuZr_Ref/'
energyFieldPrefix='pe.MD.'
strainFieldPrefix='Exy.MD.'
continuumOutPath='sct_q.out'

# Constants
kB = 8.617330350e-5
TZ = 21000
#nseed = 7536654
#np.random.seed(nseed)


"""
Reference Class
"""


class reference():

    def __init__(self, sourceDir='/Users/katianakontolati/Projects/continuum_optimization/STZ/CuZr_Ref/',
                 stressFileName='stress_rand2_matlab.txt', stressCol=6, stressDelimiter='\t', stressScaleFactor=10000,
                 maxStrain=0.5,
                 energyFieldPrefix='pe.MD.', strainFieldPrefix='Exy.MD.', nFields=101, fromContinuum=False):

        if fromContinuum:

            # Calculate T_inf
            T_ = np.fromfile('reference/tem.100', dtype=np.float32)
            ny = int(T_[0]) + 1
            nx = int(len(T_) / ny)
            T_ = T_.reshape(nx, ny)
            T_ = T_[1:, 1:]
            self.T_inf = np.amax(T_)

            self.strain = np.linspace(0, maxStrain, nFields)
            self.stress = np.zeros(nFields)

            for i in np.arange(0, nFields):
                self.stress[i] = computeContinuumDeviatoricStress('reference', i)

            self.MD_sourceDir = sourceDir
            self.MD_energyFieldPrefix = energyFieldPrefix
            self.toughness = np.trapz(self.stress, self.strain)
            self.sourceDir = 'reference/'
            self.temperatureFieldPrefix = 'tem.'
            self.strainFieldPrefix = 'Exy.'
            self.stressFieldPrefix = 'tau.'

        else:
            self.sourceDir = sourceDir
            self.energyFieldPrefix = energyFieldPrefix
            self.strainFieldPrefix = strainFieldPrefix

        self.stressDelimiter = stressDelimiter
        self.stressFileName = stressFileName
        self.stressCol = stressCol
        self.maxStrain = maxStrain
        self.stressScaleFactor = stressScaleFactor
        self.stressDelimiter = stressDelimiter
        self.nFields = nFields
        self.fromContinuum = fromContinuum

    def getStressStrainData(self):

        stressData = np.loadtxt(self.sourceDir + self.stressFileName, delimiter=self.stressDelimiter)
        self.stress = -1 * stressData[:, self.stressCol] / self.stressScaleFactor
        self.strain = np.linspace(0, self.maxStrain, len(self.stress))
        del stressData
        return self.strain, self.stress

    def getFieldValues(self, fileName, step=0, skiprows=1):
        if self.fromContinuum:
            return np.loadtxt(self.MD_sourceDir + fileName, skiprows=skiprows)
        else:
            return np.loadtxt(self.sourceDir + fileName, skiprows=skiprows)



def getContinuumField(fieldName):
    field = np.fromfile('sct_q.out' + '/' + fieldName, dtype=np.float32)
    ny = int(field[0]) + 1
    nx = int(len(field) / ny)
    field = field.reshape(nx, ny)
    field = field[1:, 1:]
    return field


def compareMatrices(A, B, stol=1e-4, method=1, \
                    scalefactor=1.0, rankreduction=True, \
                    constantrank=1, verbose=False):
    # Compute SVD of A, B:
    UA, SA, VA = np.linalg.svd(A * scalefactor, full_matrices=True)
    UB, SB, VB = np.linalg.svd(B * scalefactor, full_matrices=True)

    # Diagonalize eigenvalue arrays SA, SB:
    SA = np.diag(SA)
    SB = np.diag(SB)

    # Determine ranks of A, B:
    if constantrank == -1:
        # Rank computed to the desired eigenvalue tolerance
        rA = np.linalg.matrix_rank(SA, tol=stol)
        rB = np.linalg.matrix_rank(SB, tol=stol)
        r = np.amin([rA, rB])

    elif constantrank == 0:
        # Rank computed using smallest rank of inputs
        rA = np.linalg.matrix_rank(SA)
        rB = np.linalg.matrix_rank(SB)
        r = np.amin([rA, rB])

    elif constantrank == 1:
        # Keep the max rank of A,B
        rA = np.linalg.matrix_rank(SA, tol=stol)
        rB = np.linalg.matrix_rank(SB, tol=stol)
        r = np.amax([rA, rB])

    else:
        # Rank computed from constant rank
        # *** Need to add handling case for when requested rank exceeds input value
        if isinstance(constantrank, int):
            rA = int(constantrank)
            rB = int(constantrank)
            r = int(constantrank)
        else:
            errtxt1 = "\nWarning: Unable to deduce rank from non-integer input."
            errtxt2 = "\nUsing minimum rank instead."
            print(errtxt1 + errtxt2)

            rA = np.linalg.matrix_rank(SA)
            rB = np.linalg.matrix_rank(SB)
            r = np.amin([rA, rB])

    if r == 0:
        print("Warning: Matrix of rank 0 determined.")
        return 0, 0, 0

    # Reduce UA, UB to min or max rank:
    UA = UA[:, :r]
    UB = UB[:, :r]

    if method == 0:
        # Compute Grassmannian distance between (double infinite)  UA, UB
        R = np.dot(UA.T, UB)
        UR, SR, VR = np.linalg.svd(R, full_matrices=False)
        SR = np.round_(SR, 6)
        SR[np.where(SR > 1)] = 1.0
        theta = np.arccos(SR)
        dist = np.sqrt(abs(rA - rB) * np.pi ** 2. / 4. + \
                       np.sum(theta ** 2))
        methodtxt = "Grassmannian"

    if method == 1:
        # Compute geodesic distance on the Grassmann between U1, UB
        R = np.dot(UA.T, UB)
        UR, SR, VR = np.linalg.svd(R, full_matrices=False)
        SR = np.round_(SR, 6)
        SR[np.where(SR > 1)] = 1.0
        theta = np.arccos(SR)
        dist = np.sqrt(np.sum(theta ** 2))
        methodtxt = "Geodesic"

    # Print to screen metric and distance
    if verbose:
        text1 = "{} distance of {} determined.\n".format(methodtxt, dist)
        text2 = "Original ranks rA = {} and rB = {}\n".format(rA, rB)
        text3 = "Final rank r = {}.\n".format(r)
        print(text1 + text2 + text3)

    return dist, rA, rB


def getContinuumFieldValues(filePath, fieldName):
    """
      Retrieve continuum field quantities from binary file.
    """
    f = np.fromfile(filePath + '/' + fieldName, dtype=np.float32)
    ny = int(f[0]) + 1
    nx = int(len(f)/ny)
    f = f.reshape(nx, ny)
    f = f[1:, 1:]
    return f


def computeContinuumDeviatoricStress(filePath, frameNum):
    """
    Compute the magnitude of deviatoric shear stress for continuum output.
    """
    dev_ = getContinuumFieldValues(filePath, 'dev.{}'.format(int(frameNum)))
    dev = np.mean(dev_)

    s_ = getContinuumFieldValues(filePath, 's.{}'.format(int(frameNum)))
    s = np.mean(s_)

    tau_ = getContinuumFieldValues(filePath, 'tau.{}'.format(int(frameNum)))
    tau = np.mean(tau_)

    q_ = getContinuumFieldValues(filePath, 'q.{}'.format(int(frameNum)))
    q = np.mean(q_)

    return tau


""" 
Parameter Class
"""

class Parameter:
    name = None
    value = None
    vmin = -np.inf
    vmax = np.inf
    units_so = None
    units_hc = None
    text_so = None
    text_hc = None
    is_constant = True

    def __init__(self, name, value, vmin=-np.inf, vmax=np.inf):

        self.name = name
        self.value = value
        self.vmin = vmin
        self.vmax = vmax

        if self.vmin != -np.inf or self.vmax != np.inf:
            self.is_constant = False

        if self.name == 'beta':
            self.units_so = '1/eV'
            self.units_hc = r'$eV^{-1}$'
            self.text_so = 'beta'
            self.text_hc = r'$\beta$'

        if self.name == 'u0':
            self.units_so = 'eV'
            self.units_hc = 'eV'
            self.text_so = 'u0'
            self.text_hc = r'$u_0$'

        if self.name == 'chi_len':
            self.units_so = 'Angstroms'
            self.units_hc = r'$\AA$'
            self.text_so = 'chi_len'
            self.text_hc = r'$\ell_\chi$'

        if self.name == 'c0':
            self.units_so = '--'
            self.units_hc = '--'
            self.text_so = 'c0'
            self.text_hc = r'$c_0$'

        if self.name == 'ep':
            self.units_so = '--'
            self.units_hc = '--'
            self.text_so = 'ep'
            self.text_hc = r'$\varepsilon_0$'

        if self.name == 's_y':
            self.units_so = 'GPa'
            self.units_hc = 'GPa'
            self.text_so = 's_y'
            self.text_hc = r'$s_y$'

            try:
                assert (self.vmin >= 0)
            except:
                print("Warning! {} cannot be negative. Minimimum set to 0 K.".format(self.name))
                self.vmin = 0

            if self.name == 'chi_inf':
                self.text_so = 'chi_inf'
                self.text_hc = r'$\chi_\infty$'

            if vmin == -np.inf:
                self.vmin = 0
            elif vmin < 0:
                print("Warning! {} should not be negative. Please review value.".format(self.name))
                self.vmin = vmin
            else:
                self.vmin = vmin

            if self.vmin == -np.inf:
                self.vmin = 1
            elif self.vmin < 0:
                print("Warning! Volume cannot be less than or equal to 0.")
                print("Volume set to 1 Angstrom.")
                self.vmin = 1

    def set_min(self, value):
        self.is_constant = False
        self.vmin = value
        if self.value < self.vmin:
            self.value = self.vmin

    def set_max(self, value):
        self.is_constant = False
        self.vmax = value
        if self.value > self.vmax:
            self.value = self.vmax


"""
Parameters Class
"""


class Parameters():
    parameters = []
    ndims = None
    variableNames = []

    def __init__(self, parameters):
        self.parameters = parameters

        self.ndims = 0
        for parameter in self.parameters:
            if not parameter.is_constant:
                self.ndims += 1
                self.variableNames = np.append(self.variableNames, parameter.name)

    def list_parameters(self):
        for parameter in self.parameters:
            if parameter.is_constant:
                print("{} = {} {} and is constant.".format(
                    parameter.name, parameter.value, parameter.units_so))
            else:
                print("{} = {} {} and varies from {} to {}.".format(
                    parameter.name, parameter.value, parameter.units_so,
                    parameter.vmin, parameter.vmax))

    def getValue(self, name):
        for parameter in self.parameters:
            if parameter.name == name:
                return parameter.value

    def getVariables(self):
        output = []
        for param in self.parameters:
            if not param.is_constant:
                output = np.append(param, output)

    def setValue(self, name, value):
        for parameter in self.parameters:
            if parameter.name == name:
                parameter.value = value


def getFieldValues(fileName, step=0, skiprows=1):
    return np.loadtxt(sourceDir + fileName, skiprows=skiprows)


# Generate initial parameter list`
beta = Parameter('beta', 8.00, vmin=2, vmax=15)
u0 = Parameter('u0', -3.3626, vmin=-3.390, vmax=-3.355)
chi_len = Parameter('chi_len', 4.01, vmin=3, vmax=60)
ep = Parameter('ep', 10, vmin=5, vmax=50)
c0 = Parameter('c0', 0.3, vmin=0.05, vmax=1)
chi_inf = Parameter('chi_inf', 2730)
s_y = Parameter('s_y', 0.95, vmin=0.8, vmax=1.1)

parameters = Parameters([beta, u0, chi_inf, chi_len, c0, ep, s_y])

# Initialize dictionary to contain parameter, value pairs
initparams = {}

for param in parameters.parameters:
    initparams[param.name] = param.value

parameters.list_parameters()
print('ndims = {}'.format(parameters.ndims))

initParams = parameters

def sampling(parameters,samples,seed):
    # Sobol sequence
    initParams = parameters
    N=samples

    A = []
    for param in initParams.parameters:
        if not param.is_constant:
            A.append(Uniform(loc=param.vmin, scale=param.vmax-param.vmin))
    x1 = LHS(dist_object=A, nsamples=N, random_state=np.random.RandomState(seed), verbose=False)
    X = x1.samples
    return X

def C_matrix(A,B,k):
    C = []
    for i in range(k):
        x = np.copy(B)
        x[:,i] = A[:,i]
        C.append(x)
    return C

N = 1000  # number of samples
k = 6    # number of dimensions

# Generate samples for the two Monte Carlo simulations
A = sampling(parameters,N,700)
B = sampling(parameters,N,701)

#print(A)
#print(B)

# Generate C matrices (k in total)
C = C_matrix(A,B,k)
C1 = np.array(C[0])
C2 = np.array(C[1])
C3 = np.array(C[2])
C4 = np.array(C[3])
C5 = np.array(C[4])
C6 = np.array(C[5])


def distance(x):

    beta = x[0]
    u0 = x[1]
    chi_len = x[2]
    c0 = x[3]
    ep = x[4]
    s_y = x[5]

    if os.path.isdir('sct_q.out'):
        shutil.rmtree('sct_q.out')

    if os.path.exists('sct_q.out'):
        command = ['rm', '-r', 'sct_q.out']
        subprocess.run(command)

    PE_0 = sourceDir + energyFieldPrefix + '0'

    # Calculate final steady-state effective temperature
    PE_f = np.amax(getFieldValues('pe.MD.100'))
    T_inf = beta * (PE_f - u0) * TZ

    # Neglect b,u0 values that result in a negative (or another threshold) value of the final MD eff. temp. field
    PE_MD_final = getFieldValues('pe.MD.100')
    T_MD_final = beta * (PE_MD_final - u0) * TZ

    # Condition
    if np.any(T_MD_final < 20) == True or np.any(T_MD_final > 3000) == True:
        distance = np.nan

    command = ['shear_energy', 'qs', PE_0,
               '{}'.format(beta),
               '{}'.format(u0),
               'chi_inf', '{}'.format(T_inf),
               'chi_len', '{}'.format(chi_len),
               'c0', '{}'.format(c0),
               'ep', '{}'.format(ep),
               's_y', '{}'.format(s_y)]

    subprocess.run(command, timeout=360)

    finalCMStrainField = 'sct_q.out' + '/' + 'Exy.' + '{}'.format(100)

    if os.path.isfile(finalCMStrainField):
        strain = [0.05, 0.06, 0.08, 0.09, 0.10, 0.11, 0.12, 0.35, 0.5]
        gamOBS = np.linspace(0, 0.5, 101)

        D = []
        E = []

        for val in strain:
            gamREF = np.linspace(0, 0.5, 101)
            idREF = np.argmin(np.abs(gamREF - val))

            Aname = sourceDir + '/' + \
                    'Exy.MD.' + \
                    '{}'.format(idREF)

            D_ = np.loadtxt(Aname, skiprows=1)
            D.append(D_)

            # Get indices where strain nearest ref, obs strain
            idOBS = np.argmin(np.abs(gamOBS - val))

            Bname = 'sct_q.out' + '/' + \
                    'Exy.' + \
                    '{}'.format(idOBS)

            E_ = np.fromfile(Bname, dtype=np.float32)
            nyB = int(E_[0]) + 1
            nxB = int(len(E_) / nyB)
            E_ = E_.reshape(nxB, nyB)
            E_ = 2. * E_[1:, 1:]
            E.append(E_)

        distance_temp = np.zeros(len(D))
        r0_temp = np.zeros(len(D))
        rR_temp = np.zeros(len(D))
        for k in range(len(D)):
            distance_temp[k], r0_temp[k], rR_temp[k] = compareMatrices(D[k], E[k], stol=1e-4,
                                                                       constantrank=1, verbose=False)
        distance1 = np.mean(distance_temp)  # Grassmann distance of strain fields
        r0 = r0_temp
        rR = rR_temp

        # Also, Sum of Difference in Stress Magnitudes
        strain_cm = np.linspace(0, 0.5, 101)
        stress_cm = np.zeros(strain_cm.shape)
        distance_i = np.zeros(strain_cm.shape)

        fileName = sourceDir + '/' + 'stress_rand2_matlab.txt'
        tau_MD = np.loadtxt(fileName, delimiter="\t")
        tau_MD = -1 * tau_MD[:, 6] / 10000.
        strain_MD = np.linspace(0, 0.5, len(tau_MD))

        for i in range(0, 101):
            stress_cm[i] = s_y * computeContinuumDeviatoricStress('sct_q.out', i)

            ix = np.argmin(np.abs(strain_MD - strain_cm[i]))
            distance_i[i] = np.abs(tau_MD[ix] - stress_cm[i])

        distance_i[15:28] *= 5
        distance2 = np.sum(distance_i) / np.count_nonzero(distance_i)

        # Combined average Grassmann distance at 9 strain field snapshots and
        # average difference of stress magnitudes at all 100 frames
        distance = (1 * distance1 + 3 * distance2) / 4
        print('Distance is {}.'.format(distance))


    else:
        distance = np.nan

    return distance

# Now we have Nx(k+2) observations - model evaluations
# Each model evaluations will have (Nx1) dimensions and we will have (k+2) of those


#def Monte_Carlo(X):
#    Z = []
#    for i in range(N):
#        x = X[i, :]
#        dist = distance(x)
#        Z.append(dist)
#    return Z

#YA = Monte_Carlo(A)
#YA = np.array(YA)

#YB = Monte_Carlo(B)
#YB = np.array(YB)

#YC1 = Monte_Carlo(C1)
#YC1 = np.array(YC1)

#YC2 = Monte_Carlo(C2)
#YC2 = np.array(YC2)

#YC3 = Monte_Carlo(C3)
#YC3 = np.array(YC3)

#YC4 = Monte_Carlo(C4)
#YC4 = np.array(YC4)

#YC5 = Monte_Carlo(C5)
#YC5 = np.array(YC5)

#YC6 = Monte_Carlo(C6)
#YC6 = np.array(YC6)

YA = []
YB = []
YC1 = []
YC2 = []
YC3 = []
YC4 = []
YC5 = []
YC6 = []

for i in range(N):
    print('Iteration number:',i)
    xa = A[i, :]
    dist = distance(xa)
    YA.append(dist)

    print('Iteration number:', i)
    xb = B[i, :]
    dist = distance(xb)
    YB.append(dist)

    print('Iteration number:', i)
    xc1 = C1[i, :]
    dist = distance(xc1)
    YC1.append(dist)

    print('Iteration number:', i)
    xc2 = C2[i, :]
    dist = distance(xc2)
    YC2.append(dist)

    print('Iteration number:', i)
    xc3 = C3[i, :]
    dist = distance(xc3)
    YC3.append(dist)

    print('Iteration number:', i)
    xc4 = C4[i, :]
    dist = distance(xc4)
    YC4.append(dist)

    print('Iteration number:', i)
    xc5 = C5[i, :]
    dist = distance(xc5)
    YC5.append(dist)

    print('Iteration number:', i)
    xc6 = C6[i, :]
    dist = distance(xc6)
    YC6.append(dist)

    np.savez("Model_Output.npz", YA, YB, YC1, YC2, YC3, YC4, YC5, YC6)


# in order to open the file
# data = np.load("Model_Output.npz")
# ya = data['arr_0']   # first argument


#print(YA)
#print(YB)
#print(YC1)
#print(YC2)
#print(YC3)
#print(YC4)
#print(YC5)
#print(YC6)

