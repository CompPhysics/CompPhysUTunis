#Program which solves the one-particle Schrodinger equation 
#for a potential specified in function
#potential(). This example is for the harmonic oscillator in 3d

from  matplotlib import pyplot as plt
import numpy as np
#Function for initialization of parameters
def initialize():
    RMin = 0.0
    RMax = 10.0
    lOrbital = 0
    Dim = 600
    omega = 1.0
    return RMin, RMax, lOrbital, Dim, omega
# Here we set up the harmonic oscillator potential
def potential(r,omega):
    return r*r*omega*omega
# Here we set up the harmonic oscillator potential with interaction
def intpotential(r,omega):
    return omega*omega*r*r+1.0/r

#Get the boundary, orbital momentum and number of integration points
RMin, RMax, lOrbital, Dim, omega = initialize()

#Initialize constants
Step    = RMax/(Dim+1)
DiagConst = 2.0 / (Step*Step)
NondiagConst =  -1.0 / (Step*Step)
OrbitalFactor = lOrbital * (lOrbital + 1.0)

#Calculate array of potential values
v = np.zeros(Dim)
vint = np.zeros(Dim)
r = np.linspace(RMin,RMax,Dim)
for i in xrange(Dim):
    r[i] = RMin + (i+1) * Step;
    v[i] = potential(r[i],omega) + OrbitalFactor/(r[i]*r[i]);
    vint[i] = intpotential(r[i],omega) + OrbitalFactor/(r[i]*r[i]);

#Setting up a tridiagonal matrix and finding eigenvectors and eigenvalues
Hamiltonian = np.zeros((Dim,Dim))
IntHamiltonian = np.zeros((Dim,Dim))
Hamiltonian[0,0] = DiagConst + v[0];
Hamiltonian[0,1] = NondiagConst;
IntHamiltonian[0,0] = DiagConst + vint[0];
IntHamiltonian[0,1] = NondiagConst;
for i in xrange(1,Dim-1):
    Hamiltonian[i,i-1]  = NondiagConst;
    Hamiltonian[i,i]    = DiagConst + v[i];
    Hamiltonian[i,i+1]  = NondiagConst;
    IntHamiltonian[i,i-1]  = NondiagConst;
    IntHamiltonian[i,i]    = DiagConst + vint[i];
    IntHamiltonian[i,i+1]  = NondiagConst;

Hamiltonian[Dim-1,Dim-2] = NondiagConst;
Hamiltonian[Dim-1,Dim-1] = DiagConst + v[Dim-1];
IntHamiltonian[Dim-1,Dim-2] = NondiagConst;
IntHamiltonian[Dim-1,Dim-1] = DiagConst + vint[Dim-1];
# diagonalize and obtain eigenvalues, not necessarily sorted
EigValues, EigVectors = np.linalg.eig(Hamiltonian)
EigValuesInt, EigVectorsInt = np.linalg.eig(IntHamiltonian)
# sort eigenvectors and eigenvalues
permute = EigValues.argsort()
EigValues = EigValues[permute]
EigVectors = EigVectors[:,permute]
permute = EigValuesInt.argsort()
EigValuesInt = EigValuesInt[permute]
EigVectorsInt = EigVectorsInt[:,permute]
# now plot the results for the three lowest lying eigenstates
for i in xrange(1):
    print EigValues[i], EigValuesInt[i], 
FirstEigvector = EigVectors[:,0]
FirstEigvectorInt = EigVectorsInt[:,0]
plt.plot(r, FirstEigvector**2 ,'b-',r, FirstEigvectorInt**2 ,'g-')
plt.axis([0,3.0,0.0, 0.015])
plt.xlabel(r'$r$')
plt.ylabel(r'Radial probability $r^2|R(r)|^2$')
plt.title(r'Radial probability distributions for three lowest-lying states')
plt.savefig('eigenvectorint.pdf')
plt.show()
