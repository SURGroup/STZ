#!/usr/bin/python
import sys
import scipy.interpolate as spi

# Lists to contain all traction information
s=[];t=[]

# Lists for spline data points
es=[];et=[]

# Step size for spline data points
de=0.002

# Marker for next spline data point
ns=-1e-3

# Open file with traction information
f=open("tr_file","r");

for l in f:
    a=l.split()

    # Get strain and traction
    ss=float(a[2])
    tt=float(a[4])

    # Make main list
    s.append(ss)
    t.append(tt)

    # Make equi-sampled list
    if ss>ns:
        es.append(ss)
        et.append(tt)
        ns+=de

# Create spline
f=spi.UnivariateSpline(es,et,s=0,k=5)

# Print traction, spline-reconstructed traction, and difference
for a in range(len(s)):
    print s[a],t[a],f(s[a]),t[a]-f(s[a])
