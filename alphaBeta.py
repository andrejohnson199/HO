#!/usr/bin/env python3
import numpy as np
from flowTest import *
import matplotlib.pyplot as plt
#This testing various values for alpha and beta
#Andre Johnson

#Assign nMax,alpha, beta, gamma, ds
nMax = 6
alpha = 0.07
beta = 0.01
gamma = 1
ds = 0.01
loop = 500
nLevels = 20

#Make H,n,X
H,n = buildIntHN(nMax,alpha,beta)
X = buildX(nMax,gamma)

#Need exact H for plot
HC = HConvert(nLevels,H)
evals,evecs = np.linalg.eig(HC)
HGround = min(evals)
g = np.argmin(evals)
gvec = evecs[:,g]
XC = HConvert(nLevels,X)
numt = np.transpose(gvec)@XC@gvec
print('X=',numt)


#Flowing them
n,H,X = flowRunX(nMax,ds,loop,n,H,X)
#print(H)
print(X)

#Hold alpha, vary beta
alpha = 0
numB = 50
beta = np.linspace(0,0.1,numB)
H00 = np.zeros(numB)
HGrounds = np.zeros(numB)
AllH = np.zeros((numB,nLevels))
beta2 = -6*beta**2
for i in range(numB):
    #Need to make H,n 
    H,n = buildIntHN(nMax,alpha,beta[i])
    #Get exact H
    HC = HConvert(nLevels,H)
    evals,evecs = np.linalg.eigh(HC)
    HGrounds[i] = min(evals)
    AllH[i,:] = evals
    #Flow
    n, H = flowRun(nMax,ds,loop,n,H)
    H00[i] = H[0][0]

#Plotting
legends = ['Flow','Exact','Pert']
plt.plot(beta,H00)
plt.plot(beta,HGrounds)
plt.plot(beta,beta2)
plt.xlabel('beta')
plt.ylabel('Ground State')
plt.legend(legends)
plt.savefig('betaRun.pdf')

#Plotting All eigenvalues vs beta
plt.figure(2)
plt.plot(beta,AllH)
plt.xlabel('beta')
plt.ylabel('Energies')
plt.savefig('energyRunBeta.pdf')

#Hold Beta run alpha
beta = 0.01
numA = 50
alpha = np.linspace(0,0.2,numA)
H00 = np.zeros(numA)
HGrounds = np.zeros(numA)
AllH = np.zeros((numA,nLevels))
pert = -6*beta**2-2*alpha**2
for i in range(numA):
    #Need to make H,n 
    H,n = buildIntHN(nMax,alpha[i],beta)
    #Get exact H
    HC = HConvert(nLevels,H)
    evals,evecs = np.linalg.eigh(HC)
    HGrounds[i] = min(evals)
    AllH[i,:] = evals

    #Flow
    n, H = flowRun(nMax,ds,loop,n,H)
    H00[i] = H[0][0]

#Plotting
plt.figure(3)
legends = ['Flow','Exact','pert']
plt.plot(alpha,H00)
plt.plot(alpha,HGrounds)
plt.plot(alpha,pert)
plt.xlabel('alpha')
plt.ylabel('Ground State')
plt.legend(legends)
plt.savefig('alphaRun.pdf')

#Plotting all energies for alpha
plt.figure(4)
plt.plot(alpha,AllH)
plt.xlabel('alpha')
plt.ylabel('Energies')
plt.savefig('energyRunAlpha.pdf')

#Testing 1d wavefunctions


