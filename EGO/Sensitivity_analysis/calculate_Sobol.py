#!/usr/bin/env python

import numpy as np

# Load data produced by 'Sobol_indices.py'
data = np.load("Model_Output.npz")
ya = data['arr_0']   # first argument
yb = data['arr_1']
yc1 = data['arr_2']
yc2 = data['arr_3']
yc3 = data['arr_4']
yc4 = data['arr_5']
yc5 = data['arr_6']
yc6 = data['arr_7']

# Save them in a single list
Y = [ya,yb,yc1,yc2,yc3,yc4,yc5,yc6]

# Assign np.inf values to all the np.nan values
for i in range(len(Y)):
    Y[i][np.isnan(Y[i])] = np.inf
    #print('for the following indices:',np.argwhere(np.isnan(Y[i])).T)

# Save in separate lists
ya_ = Y[0]
yb_ = Y[1]
yc1_ = Y[2]
yc2_ = Y[3]
yc3_ = Y[4]
yc4_ = Y[5]
yc5_ = Y[6]
yc6_ = Y[7]

# Keep the matrix indices where all 8 matrices have in common non-inf values
n = len(ya_)
indices = []
counter = 0
for i in range(n):
    if ya_[i] != np.inf and yb_[i] != np.inf and yc1_[i] != np.inf and yc2_[i] != np.inf and yc3_[i] != np.inf and yc4_[i] != np.inf and yc5_[i] != np.inf and yc6_[i] != np.inf:
        counter = counter + 1
        indices.append(i)

# make indices an numpy array
indices = np.array((indices))
print('The number of non-nans is: {} out of the {}'.format(counter, len(ya_)))
#print(indices)

# Produce final Monte Carlo output
ya_f = ya_[indices]
yb_f = yb_[indices]
yc1_f = yc1_[indices]
yc2_f = yc2_[indices]
yc3_f = yc3_[indices]
yc4_f = yc4_[indices]
yc5_f = yc5_[indices]
yc6_f = yc6_[indices]

Cf = [yc1_f,yc2_f,yc3_f,yc4_f,yc5_f,yc6_f]

# Compute
# total number of observations to be considered for the calculation of first-order sensitivity indices
N = len(ya_f)
# mean value
f_02 = ((ya_f.sum())/N)**2

den = (1/N)*((ya_f*ya_f).sum()) - f_02  # same denominator for all indices
print('denominator is:',den)

# Compute first-order sensitivity indices
total = 0
for i in range(len(Cf)):
    c = np.array((Cf[i]))
    num = (1/N)*((ya_f*c).sum()) - f_02
    S = num/den
    print('First-order Sobol index no. {} is equal to {}'.format((i+1),S))
    total = total + S

print('sum of indices is:',total)

# Compute total-effect sensitivity indices
total_t = 0
for i in range(len(Cf)):
    c = np.array((Cf[i]))
    num = (1/N)*((yb_f*c).sum()) - f_02
    S_T = 1 - (num/den)
    print('Total-effect Sobol index no. {} is equal to {}'.format((i+1),S_T))
    total_t = total + S

print('sum of total indices is:',total_t)
