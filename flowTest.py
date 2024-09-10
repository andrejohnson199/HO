#!/usr/bin/env python3
import math
import numpy as np
#Test run for all possible n and H combinations
#Andre Johnson

#N Body max
#nMax = 4
def flowCombo(nMax,n,H):
    dH = np.zeros((nMax+1,nMax+1))
    #Outer loop for nBody
    for nBody in range(nMax+1):
        #print('Body Level: ',nBody)
        #Number of contractions between n and H
        for nContract in range(1,nMax+1):
            #print(nContract)
            # n index
            for x in range(nMax-nContract+1):
                a = nContract + x
                # H first index
                for y in range(nMax-nContract+1):
                    b = nContract + y
                    # H second index
                    c = nBody - x - y
                    if c > -1 and b+c <= nMax:
                        #Finding dh element
                        e1 = x + c
                        e2 = y
                        if e2 > e1:
                            temp = e1
                            e1 = e2
                            e2 = temp
                        #if e1 == e2: #This is for the "even" diagrams that allow n and H to flip: 00,11,22  
                         #   print('n:',a,', H:',b,',',c,', dH element:',e1,e2,', X2')
                        #else:
                         #   print('n:',a,', H:',b,',',c,', dH element:',e1,e2)

                        #Adding element to dH
                        flowValue = update(nContract,x,y)*n[a]*H[b][c]
                        if e1 == e2: flowValue+=flowValue
                        dH[e1][e2] += flowValue
                        dH[e2][e1] = dH[e1][e2]

                        #FIXING DH 40

                        #Check on the H00 elements
                        #if e1==0 and e2==0:
                            #print(update(nContract,x,y),n[a],H[b][c])
    
    return dH

#Giving the resulting dh element value
def update(i,x,y):
    return -1*math.factorial(i)*math.comb(i+x,x)*math.comb(i+y,y)
#testing = flowCombo(4,0,0)

#Build the x^3 + x^4 HO
def buildIntH(nMax,alpha,beta):
    H = np.zeros((nMax+1,nMax+1))
    H[1][1] = 1
    for col in range(1,nMax+1):
        for row in range(col+1):
            if row+col == 3:
                H[row][col] = alpha * math.comb(row+col,col)
                H[col][row] = H[row][col]
            if row+col == 4:
                H[row][col] = beta * math.comb(row+col,col)
                H[col][row] = H[row][col]
    return H

#Same as buildIntH but also returns n
def buildIntHN(nMax,alpha,beta):
    H = buildIntH(nMax,alpha,beta)
    n = setN(nMax,H)
    return H,n

#Builds the x coeff matrix
def buildX(nMax,gamma):
    X = np.zeros((nMax+1,nMax+1))
    X[0][1] = 1*gamma
    X[1][0] = X[0][1]
    return X

#Loop flow method
#Will switch over to while loop eventually
def flowRunX(nMax,ds,loop,n,H,X):
    newH = np.copy(H)
    newN = np.copy(n)
    newX = np.copy(X)
    #Looping
    for i in range(loop):
        dH = flowCombo(nMax,newN,newH)
        dX = flowCombo(nMax,newN,newX)
        newH += dH*ds
        newX += dX*ds
        newN = setN(nMax,newH)
    return newN,newH,newX

#Loop flow method
#Will switch over to while loop eventually
def flowRun(nMax,ds,loop,n,H):
    newH = np.copy(H)
    newN = np.copy(n)
    #Looping
    for i in range(loop):
        dH = flowCombo(nMax,newN,newH)
        newH += dH*ds
        newN = setN(nMax,newH)
    return newN,newH

#Set n method 
def setN(nMax,H):
    n = np.zeros(nMax+1)
    for i in range(nMax+1):
        n[i] = H[0][i]
    return n

#Convert H coeff matrix into Hamiltonian Matrix
def HConvert(nLevels,hMat):
    reMat = np.zeros((nLevels,nLevels))
    rows,columns = hMat.shape
    for i in range(rows):
        for j in range(columns):
            reMat+=hMat[i][j]*HOBuilder(nLevels,i,j)
    return reMat

#Method for building a delta matrix
def delta(nLevels,x,y):
    d = np.zeros((nLevels,nLevels))
    for i in range(nLevels):
        for j in range(nLevels):
            if i+x == j+y:
                d[i][j] = 1
    return d

#Method for building creation/annihilation in HO
def HOBuilder(nLevels,nAt,nA):
    delMat = delta(nLevels,0,nAt-nA)
    return HOdelta(nLevels,delMat,nAt,nA)

#Adds coeff to delta matrix
def HOdelta(nLevels,dMat,nAt,nA):
    for i in range(nLevels):
        for j in range(nLevels):
            if dMat[i][j] != 0:
                dMat[i][j] *= coeffBuild(j,-nA)*coeffBuild(j-nA,nAt)
    return dMat

#Returns coeff from creation/annihilation
def coeffBuild(ket,n):
    #Ket could be negative
    if ket < 0:
        return 0
    #empty
    if n==0: return 1
    returnVal = 1
    #Annihilation
    if n < 0:
        for i in range(-n):
            if ket-i<0:
                returnVal = 0
            else:
                returnVal *= (ket-i)**(1/2)
    else:
        for i in range(1,n+1):
            returnVal *= (ket+i)**(1/2)
    return returnVal

    
