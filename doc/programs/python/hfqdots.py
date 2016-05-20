import numpy as np
from decimal import Decimal
# expectation value for the one body part, Harmonic oscillator in two dimensions
def onebody(i, n, m):
        homega = 1.0
        return homega*(2*n[i] + m[i] + 1)

if __name__ == '__main__':

        Nparticles = 2
        """ Read quantum numbers from file """
        index = []
        n = []
        mj = []
        spOrbitals = 0
        with open("qdnumbers.dat", "r") as qnumfile:
                for line in qnumfile:
                        nums = line.split()
                        if len(nums) != 0:
                                index.append(int(nums[0]))
                                n.append(int(nums[1]))
                                mj.append(int(nums[4]))
                                spOrbitals += 1


        """ Read two-electron interaction elements (integrals) from file, brute force 4-dim array """
        interaction = np.zeros([spOrbitals, spOrbitals, spOrbitals, spOrbitals])
        with open("twobody.dat", "r") as infile:
                for line in infile:
                        number = line.split()
                        a = int(number[0]) - 1
        b = int(number[1]) - 1
                        c = int(number[2]) - 1
                        d = int(number[3]) - 1
                        #print a, b, c, d, float(l[4])
                        interaction[a][b][c][d] = Decimal(number[4])
        """ Set up single-particle integral """
        singleparticleH = np.zeros(spOrbitals)
       for i in range(spOrbitals):
                singleparticleH[i] = Decimal(onebody(i, n, m))

        """ Star HF-iterations, preparing variables and density matrix """

        """ Coefficients for setting up density matrix, assuming only one along the diagonals """
        C = np.eye(spOrbitals) # HF coefficients
        DensityMatrix = np.zeros([spOrbitals,spOrbitals])
        for gamma in range(spOrbitals):
            for delta in range(spOrbitals):
                sum = 0.0
                for i in range(Nparticles):
                    sum += C[gamma][i]*C[delta][i]
                DensityMatrix[gamma][delta] = Decimal(sum)
        maxHFiter = 100
        epsilon =  1.0e-5
        difference = 1.0
        hf_count = 0
        oldenergies = np.zeros(spOrbitals)
        newenergies = np.zeros(spOrbitals)
        while hf_count < maxHFiter and difference > epsilon:
                print "############### Iteration %i ###############" % hf_count
                HFmatrix = np.zeros([spOrbitals,spOrbitals])
                for alpha in range(spOrbitals):
                        for beta in range(spOrbitals):
                            """  If tests for two-dimensional system """
                            if mj[alpha] != mj[beta]: continue
                            """  Setting up the Fock matrix using the density matrix and antisymmetrized  in\
teraction in m-scheme """
                            sumFockTerm = 0.0
                            for gamma in range(spOrbitals):
                                for delta in range(spOrbitals):
                                    if (mj[alpha]+mj[gamma]) != (mj[beta]+mj[delta]): continue
                                    sumFockTerm += DensityMatrix[gamma][delta]*interaction[alpha][gamma][b\
eta][delta]
                            HFmatrix[alpha][beta] = Decimal(sumFockTerm)
                            """  Adding the one-body term, here plain harmonic oscillator """
                            if beta == alpha:   HFmatrix[alpha][alpha] += singleparticleH[alpha]
                spenergies, C = np.linalg.eigh(HFmatrix)
                """ Setting up new density matrix in m-scheme """
                DensityMatrix = np.zeros([spOrbitals,spOrbitals])
                for gamma in range(spOrbitals):
                    for delta in range(spOrbitals):
                        sum = 0.0
                        for i in range(Nparticles):
                            sum += C[gamma][i]*C[delta][i]
       DensityMatrix[gamma][delta] = Decimal(sum)
                newenergies = spenergies
""" Brute force computation of difference between previous and new sp HF energies """
                sum =0.0
                for i in range(spOrbitals):
                    sum += (abs(newenergies[i]-oldenergies[i]))/spOrbitals
                difference = sum
                oldenergies = newenergies
                print "Single-particle energies, ordering may have changed "
                for i in range(spOrbitals):
                    print('{0:4d}  {1:.4f}'.format(i, Decimal(oldenergies[i])))
                hf_count += 1
