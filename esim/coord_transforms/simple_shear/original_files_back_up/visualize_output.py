import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
#plt.switch_backend('Qt5Agg')
#plt.ion()
import os, sys, fnmatch, subprocess

#convert binary array to (nx times ny) matrix
def Snapshot_to_field(filename):
	# Final CM plastic strain
	S = np.fromfile('sct_q.out/'+filename, dtype=np.float32)
	ny = int(S[0]) + 1
	nx = int(len(S)/ny)
	S = S.reshape(nx, ny)
	S = S[1:, 1:] #remove first row and column
	return S

#os.mkdir('figures_Exy')
#os.mkdir('figures_tem')
#os.mkdir('matrices')

T_mean = np.arange(200,400,100)
for T in T_mean:
    command = ['./shear_no_MD', 'T_mean', '{}'.format(T)]
    subprocess.run(command, shell=True)
    os.makedirs('T_mean/{}'.format(T))
    
    # Outputs
    s = np.arange(0, 101)
    field1 = 'Exy'
    field2 = 'tem'
    
    for s_ in s:
        filename = '{}.{}'.format(field1, s_)
        F = Snapshot_to_field(filename)
        filename = '{}.{}'.format(field2, s_)
        P = Snapshot_to_field(filename)

        # Shear stress Exy
        fig1 = plt.figure()
        plt.contourf(F)
        plt.colorbar()
        plt.savefig('T_mean/{}/{}.{}.png'.format(T,field1, s_))
        plt.close()
    	
        # Temperature
        fig2 = plt.figure()
        plt.contourf(P)
        plt.colorbar()
        plt.savefig('T_mean/{}/{}.{}.png'.format(T,field2, s_))
        plt.close()
    
        np.savez('T_mean/{}/{}.{}.npz'.format(T,field1, s_), F=F)
   
