#!/home/flavio/anaconda3/bin/python
from __future__ import print_function
import math
import random
import sys
import numpy as np
from scipy import linalg as la

__author__ = "Flavio Riche"
__version__ = "1.0"

# -----------------------------------------------------------------------------------------
# This script models the AA model. It diagonalizes the Hamiltonian for a system of finite size L,
# and twisted boundary conditions. The script computes and prints the ipr and the energy spectrum
# for a list of potential values ranging from 0 to 4.
# -----------------------------------------------------------------------------------------


#parameters

L=100 #System size
tau= (5 ** 0.5 - 1) / 2 #Inverse of the Golden Ratio
k = np.random.ranf() * 2 * np.pi #parameter for twisted boundary conditions
phi = np.random.ranf() * 2 * np.pi #random phase for the AA potential 

# Generate the Hamiltonian matrix 
def generate_hamiltonian(L,V):
    H = np.zeros((L, L),dtype='complex_')
    for i in range(L):
        ip1 = (i + 1) % L
# on-site quasi-periodic potential
        H[i, i] = V*np.cos(2. * np.pi * tau * i + phi)
# first-neighbor hopping -- twisted boundary conditions
        H[i, ip1] = np.exp(1j * k)
        H[ip1, i] = np.exp(-1j * k)
    return H

# Calculate the IPR
def compute_ipr(eigenstates):
  ipr=np.sum(abs(eigenstates)**4,axis=0)
  return ipr

# Diagonalization
V_list=np.arange(0,4.1,0.1)
for V in V_list:
    H=generate_hamiltonian(L,V)
    (energy_levels,eigenstates)=la.eigh(H)
    idx = energy_levels.argsort()[::1] 
    energy_levels=energy_levels[idx]
    eigenstates=eigenstates[:,idx]
# Calculation of IPR
    ipr=compute_ipr(eigenstates) 
# print V, IPR and energy spectrum
    #filename='energy_spectrum_.dat'
    #f=open(filename,'w')
    #for j in range(L):
        print(str(V)+" "+str(ipr[j])+" "+str(energy_levels[j])+"\n")
    #    f.write(str(V)+" "+str(ipr[j])+" "+str(energy_levels[j])+"\n")
    #f.close()
