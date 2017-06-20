import os
import sys
import shutil
import numpy as np
import random
from time import time
import multiprocessing
import datetime
import math
from scipy.stats import norm
import matplotlib.pyplot as plt
from collections import defaultdict
from operator import itemgetter
import csv
from pdbReader import PDBReader
from pdbAligner import PDBAligner

def evaluator(x,refPosAtoms, modPosAtoms):
    translation = []
    rotation = []

    translation = [ 200*i-100 for i in x[:3] ]
    rotation = [ 2*math.pi*i for i in x[3:] ]

    aligner = PDBAligner()
    solution = aligner.transform( modPosAtoms, translation, rotation )

    f = aligner.calcRMSD( refPosAtoms, solution )
    res = {'x1':f-5,'x2':2*f, 'x3':f-5,'x4':2*f, 'x5':f-5,'x6':2*f }
    fitness = dict(Obj=f,**res)
    return fitness

def mp_evaluator(x, refPosAtoms, modPosAtoms):
    nprocs = 4
    # create pool
    pool = multiprocessing.Pool(processes=nprocs)
    results = [pool.apply_async(evaluator,[c, refPosAtoms, modPosAtoms]) for c in x]
    pool.close()
    pool.join()
    #print "results: ", results
    f = [r.get()['Obj'] for r in results]
    #print f
    for r in results:
        del r.get()['Obj']
    # maximization or minimization problem
    maximize = False
    return (f, [r.get() for r in results],maximize)

def initialize(ants,var):
    X = np.random.uniform(low=0,high=1,size=(ants,var))
    return X

def evolve( refPosAtoms, modPosAtoms):
    start_time = time()

    # number of variables
    parameters_v = ['tX', 'tY', 'tZ', 'rX', 'rY', 'rZ']
    response_v = ['x1','x2', 'x3', 'x4', 'x5', 'x6']

    # number of variables
    numVariables = len(parameters_v)
    # size of solution archive
    solutionSize = 100
    # number of ants
    nAnts = 50

    # parameter q
    q = 0.02

    # standard deviation
    qk = q*solutionSize

    # parameter xi (like pheromone evaporation)
    xi = 0.8

    # maximum iterations
    maxiter = 10

    # bounds of variables
    Up = [1]*numVariables
    Lo = [0]*numVariables

    # initilize matrices
    S = np.zeros((solutionSize,numVariables))
    S_f = np.zeros((solutionSize,1))

    # inicializar matriz de solucoes
    Srand = initialize(solutionSize,numVariables)

    f,S_r,maximize = mp_evaluator(Srand, refPosAtoms, modPosAtoms)

    S_responses = []

    for i in range(len(S_r)):
        S_f[i] = f[i]
        k = S_r[i]
        row = []
        for r in response_v:
            row.append( k[r])
        S_responses.append(row)

    # add responses and "fitness" column to solution
    S = np.hstack((Srand,S_responses,S_f))
    # sort according to fitness (last column)
    S = sorted(S, key=lambda row: row[-1],reverse = maximize)
    S = np.array(S)

    # initilize weight array with pdf function
    w = np.zeros((solutionSize))
    for i in range(solutionSize):
        w[i] = 1/(qk*math.sqrt(2*math.pi))*math.exp(-math.pow(i,2)/(2*math.pow(q,2)*math.pow(solutionSize,2)))

    # initialize variables
    iterations = 1
    best_par = []
    best_obj = []
    best_sol = []
    best_res = []
    worst_obj = []
    best_par.append(S[0][:numVariables])
    best_obj.append(S[0][-1])
    best_sol.append(S[0][:])
    best_res.append(S[0][numVariables:-1])
    worst_obj.append(S[-1][-1])

    stop = 0

    # iterations
    while True:
        # choose Gaussian function to compose Gaussian kernel
        p = w/sum(w)

        # find best and index of best
        max_prospect = np.amax(p)
        ix_prospect = np.argmax(p)
        selection = ix_prospect

        # calculation of G_i
        # find standard deviation sigma
        sigma_s = np.zeros((numVariables,1))
        sigma = np.zeros((numVariables,1))
        for i in range(numVariables):
            for j in range(solutionSize):
                sigma_s[i] = sigma_s[i] + abs(S[j][i] - S[selection][i])
            sigma[i] = xi / (solutionSize -1) * sigma_s[i]

        Stemp = np.zeros((nAnts,numVariables))
        ffeval = np.zeros((nAnts,1))
        res = np.zeros((nAnts,len(response_v)))
        for k in range(nAnts):
            for i in range(numVariables):
                Stemp[k][i] = sigma[i] * np.random.random_sample() + S[selection][i]
                if Stemp[k][i] > Up[i]:
                    Stemp[k][i] = Up[i]
                elif Stemp[k][i] < Lo[i]:
                    Stemp[k][i] = Lo[i]
        #print Stemp
        f,S_r,maximize = mp_evaluator(Stemp, refPosAtoms, modPosAtoms)

        S_f = np.zeros((nAnts,1))
        S_responses = []

        for i in range(len(S_r)):
            S_f[i] = f[i]
            k = S_r[i]
            row = []
            for r in response_v:
                row.append(k[r])
            S_responses.append(row)

        # add responses and "fitness" column to solution
        Ssample = np.hstack((Stemp,S_responses,S_f))

        # add new solutions in the solutions table
        Solution_temp = np.vstack((S,Ssample))

        # sort according to "fitness"
        Solution_temp = sorted(Solution_temp, key=lambda row: row[-1],reverse = maximize)
        Solution_temp = np.array(Solution_temp)

        #print Solution_temp
        # keep best solutions
        S = Solution_temp[:solutionSize][:]

        # keep best after each iteration
        best_par.append(S[0][:numVariables])
        best_obj.append(S[0][-1])
        best_res.append(S[0][numVariables:-1])
        best_sol.append(S[0][:])
        worst_obj.append(S[-1][-1])

        print best_sol[0][-1]
        iterations += 1
        if iterations > maxiter or stop > 5 or best_sol[0][-1] == 0.0:
            break

    total_time_s = time() - start_time
    total_time = datetime.timedelta(seconds=total_time_s)
    #total_time = formatTD(total_time)

    best_sol = sorted(best_sol, key=lambda row: row[-1],reverse = maximize)

    print best_sol[0][-1]

def main():
    refPdb = PDBReader( "files/reference.pdb" )
    modPdb = PDBReader( "files/reference.pdb" )

    oldRef = refPdb.posAtoms
    oldMod = modPdb.posAtoms

    modPdb.adjustAtoms( refPdb.atoms, refPdb.aminoAcids )

    refPdb.calcBackbonePos()
    modPdb.calcBackbonePos()
    refPdb.calcCaPos()
    modPdb.calcCaPos()

    '''print "#########All atoms#########"
    for i in range( 30 ):
        evolve( refPdb.posAtoms, modPdb.posAtoms )

    print "#########Backbone atoms#########"
    for i in range( 30 ):
        evolve( refPdb.backbone, modPdb.backbone )'''

    print "#########Alpha atoms#########"
    for i in range( 1 ):
        evolve( refPdb.alpha, modPdb.alpha )

if __name__ == "__main__":
    main()