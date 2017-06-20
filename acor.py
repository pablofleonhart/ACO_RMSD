'''
    ==============================================================
    Ant Colony Optimization algorithm for continuous domains ACO_R
    ==============================================================
    author: Andreas Tsichritzis <tsadreas@gmail.com>
'''

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
import copy
import csv
from pdbReader import PDBReader
from pdbAligner import PDBAligner
from testAligner import tAligner

path = os.getcwd()
# get folder name from input
case = str(sys.argv[1])
if os.path.exists(path+'/'+case):
    raw_input('Folder exists! Press enter to delete it and continue')
    shutil.rmtree(path+'/'+case)
os.mkdir(path+'/'+case)
os.chdir(path+'/'+case)

def evaluator(x,refPosAtoms, modPosAtoms):
    '''Evaluator function, returns fitness and responses values'''
    # give the normalized candidates values inside the real design space
    translation = []
    rotation = []

    translation = [ 200*i-100 for i in x[:3] ]
    rotation = [ 2*math.pi*i for i in x[3:] ]

    transformation = [translation[0],translation[1],translation[2],rotation[0],rotation[1],rotation[2]]
    #print "tr", translation
    #print "rt", rotation
    #print transformation
    aligner = PDBAligner()
    solution = aligner.transform( copy.deepcopy( modPosAtoms ), translation, rotation )

    f = aligner.calcRMSD( refPosAtoms, solution )
    #print "f", f

    ta = tAligner()
    f = ta.rmsd( transformation, refPosAtoms, copy.deepcopy( modPosAtoms ) )
    # calculate values for other responses
    res = {'x1':f,'x2':f, 'x3':f,'x4':f, 'x5':f,'x6':f }
    fitness = dict(Obj=f,**res)
    return fitness

def mp_evaluator(x,refPosAtoms, modPosAtoms):
    '''Multiprocessing evaluation'''
    # ste number of cpus
    nprocs = 1
    # create pool
    pool = multiprocessing.Pool(processes=nprocs)
    results = [pool.apply_async(evaluator,[c,refPosAtoms, modPosAtoms]) for c in x]
    pool.close()
    pool.join()
    f = [r.get()['Obj'] for r in results]
    for r in results:
        del r.get()['Obj']
    # maximization or minimization problem
    maximize = False
    return (f, [r.get() for r in results],maximize)

def initialize(ants,var):
    '''Create initial solution matrix'''
    X = np.random.uniform(low=0,high=1,size=(ants,var))
    '''for i in xrange( ants ):
        X[i][0] = 0.5;
        X[i][1] = 0.5;
        X[i][2] = 0.5;
        X[i][3] = 0;
        X[i][4] = 0;
        X[i][5] = 0;
    print "X", X'''
    return X

def init_observer(filename,matrix,parameters,responses):
    '''Initial population observer'''
    p = []
    r = []
    f = []
    res = ['{0:>10}'.format(i)[:10] for i in responses]
    par = ['{0:>10}'.format(i)[:10] for i in parameters]
    for i in range(len(matrix)):
        p.append(matrix[i][0:len(parameters)])
        r.append(matrix[i][len(parameters):-1])
        f.append(matrix[i][-1])
    r = np.array(r)
    p = np.array(p)

    for i in range(len(r)):
        r[i] = ['{0:>10}'.format(r[i][j])[:10] for j in range(len(responses))]

    for i in range(len(p)):
        p[i] = ['{0:>10}'.format(p[i][j])[:10] for j in range(len(parameters))]

    f = ['{0:>10}'.format(i)[:10] for i in f]

    iteration = 0

    filename.write('{0:>10}, {1}, {2:>10}, {3}\n'.format('Iteration',', '.join(map(str, par)),'Fitness',', '.join(map(str, res))))

    for i in range(len(matrix)):
        filename.write('{0:>10}, {1}, {2:>10}, {3}\n'.format(iteration,', '.join(map(str, p[i])),f[i],', '.join(map(str, r[i]))))

def iter_observer(filename,matrix,parameters,responses,iteration):
    '''Iterations observer'''
    p = []
    r = []
    f = []
    for i in range(len(matrix)):
        p.append(matrix[i][0:len(parameters)])
        r.append(matrix[i][len(parameters):-1])
        f.append(matrix[i][-1])
    r = np.array(r)
    p = np.array(p)

    for i in range(len(r)):
        r[i] = ['{0:>10}'.format(r[i][j])[:10] for j in range(len(responses))]

    for i in range(len(p)):
        p[i] = ['{0:>10}'.format(p[i][j])[:10] for j in range(len(parameters))]

    f = ['{0:>10}'.format(i)[:10] for i in f]

    for i in range(len(matrix)):
        filename.write('{0:>10}, {1}, {2:>10}, {3}\n'.format(iteration,', '.join(map(str, p[i])),f[i],', '.join(map(str, r[i]))))

def correct_par(filename,par):
    """Replace normalized values with real"""
    columns = defaultdict(list)
    with open(filename) as f:
        reader = csv.DictReader(f,skipinitialspace=True)
        for row in reader:
            for (k,v) in row.items():
                columns[k].append(v)
        keys = columns.keys()
        for p in par:
            if p in keys:
                col = []
                for i,k in enumerate(columns[p]):
                    k = float(k)
                    if p in par:
                        n = 10*k-5
                    col.append(n)
                columns[p] = col

    outputfile = filename

    file = open(outputfile,'w+')
    head = []
    head.append('Iteration')
    for i in par:
        head.append(i)
    head.append('Fitness')
    for i in keys:
        if i not in head:
            head.append(i)
    par = ['{0:>10}'.format(i)[:10] for i in par]
    line = ['{0:>10}'.format(l)[:10] for l in head]
    file.write('{0}\n'.format(', '.join(map(str, line))))
    for i in range(len(columns.get('Iteration'))):
        line = []
        for j in head:
            line.append(columns.get(j)[i])
        line = ['{0:>10}'.format(l)[:10] for l in line]
        file.write('{0}\n'.format(', '.join(map(str, line))))
    file.close()

def formatTD(td):
    """ Format time output for report"""
    days = td.days
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%s days %s h %s m %s s' % (days, hours, minutes, seconds)

def evolve(refPosAtoms, modPosAtoms,display):
    '''Executes the optimization'''
    start_time = time()

    # number of variables
    parameters_v = ['tX', 'tY', 'tZ', 'rX', 'rY', 'rZ']
    response_v = ['x1','x2', 'x3', 'x4', 'x5', 'x6']

    # create output file
    projdir = os.getcwd()
    ind_file_name = '{0}/results.csv'.format(projdir)
    ind_file = open(ind_file_name, 'w')

    # number of variables
    nVar = len(parameters_v)
    # size of solution archive
    nSize = 200
    # number of ants
    nAnts = 150

    # parameter q
    q = 0.0001

    # standard deviation
    qk = q*nSize

    # parameter xi (like pheromone evaporation)
    xi = 0.85

    # maximum iterations
    maxiter = 2000
    # tolerance
    errormin = 0.01

    # bounds of variables
    Up = [1]*nVar
    Lo = [0]*nVar

    # initilize matrices
    S = np.zeros((nSize,nVar))
    S_f = np.zeros((nSize,1))

    plt.figure()

    # initialize the solution table with uniform random distribution and sort it
    print '-----------------------------------------'
    print 'Starting initilization of solution matrix'
    print '-----------------------------------------'

    Srand = initialize(nSize,nVar)
    print "Inicial", Srand
    f,S_r,maximize = mp_evaluator(Srand, refPosAtoms, modPosAtoms)
    print f

    S_responses = []

    for i in range(len(S_r)):
        S_f[i] = f[i]
        k = S_r[i]
        row = []
        for r in response_v:
            row.append(k[r])
        S_responses.append(row)

    # add responses and "fitness" column to solution
    S = np.hstack((Srand,S_responses,S_f))
    # sort according to fitness (last column)
    S = sorted(S, key=lambda row: row[-1],reverse = maximize)
    S = np.array(S)
    init_observer(ind_file,S,parameters_v,response_v)

    # initilize weight array with pdf function
    w = np.zeros((nSize))
    for i in range(nSize):
        w[i] = ( 1/( qk * math.sqrt( 2*math.pi ) ) )*math.exp( -math.pow( i, 2 )/( 2*math.pow( q, 2 )*math.pow( nSize, 2 ) ) )
        #print i, w[i]

    if display:
        x = []
        y = []
        for i in S:
            x.append(i[0])
            y.append(i[1])

        plt.scatter(x,y)
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.pause(2)
        plt.cla()

    # initialize variables
    iterations = 1
    best_par = []
    best_obj = []
    best_sol = []
    best_res = []
    worst_obj = []
    best_par.append(S[0][:nVar])
    best_obj.append(S[0][-1])
    best_sol.append(S[0][:])
    best_res.append(S[0][nVar:-1])
    worst_obj.append(S[-1][-1])

    p = w/sum(w)
    stop = 0

    # iterations
    while True:
        print '-----------------------------------------'
        print 'Iteration', iterations
        print '-----------------------------------------'
        # choose Gaussian function to compose Gaussian kernel        
        #print "w", w
        #print "p", p
        # find best and index of best
        '''max_prospect = np.amax(p)
        ix_prospect = np.argmax(p)
        selection = ix_prospect
        print "selection", selection
        for k in range( nAnts ):
            # escolhe uma solucao
            cs = np.random.random_sample()
            #print cs
            total = 0
            for z in xrange( nSize-1, -1, -1 ):
                total += p[z]
                if cs <= total:
                    sol = z
                    break
        # calculation of G_i
        # find standard deviation sigma
        sigma_s = np.zeros((nVar,1))
        sigma = np.zeros((nVar,1))
        #print "sigma", sigma
        for i in range(nVar):
            for j in range(nSize):
                sigma_s[i] += abs(S[j][i] - S[sol][i])
            sigma[i] = xi * ( sigma_s[i] / (nSize -1) )
        Stemp = np.zeros((nAnts,nVar))
        ffeval = np.zeros((nAnts,1))
        res = np.zeros((nAnts,len(response_v)))
        '''
        Stemp = np.zeros((nAnts,nVar))

        # para cada formiga
        for k in range( nAnts ):
            # escolhe uma solucao
            cs = np.random.random_sample()
            #print cs
            total = 0
            for z in xrange( nSize-1, -1, -1 ):
                total += p[z]
                if cs <= total:
                    sol = z
                    break

            #print sol

            # para cada variavel
            for i in range( nVar ):
                # calcula desvio padrao da solucao 'sol'
                sigma = 0
                for y in xrange( nSize ):
                    sigma += abs( S[y][i] - S[sol][i] )/( nSize-1 )

                #print "sigma", i, "=", sigma

                # calcula valor de 'i' com as funcoes gaussianas
                x = np.random.random_sample()
                gi = w[sol]*math.exp( -math.pow( x - S[sol][i], 2 ) / ( 2*math.pow( sigma, 2 ) ) )* (1/( sigma*math.pow( 2*math.pi, 2 ) ))

                #print gi
                #Stemp[k][i] = sigma[i] * x + S[selection][i]
                Stemp[k][i] = gi
                if Stemp[k][i] > Up[i]:
                    Stemp[k][i] = Up[i]
                elif Stemp[k][i] < Lo[i]:
                    Stemp[k][i] = Lo[i]

        # find best and index of best
        '''
        max_prospect = np.amax(p)
        ix_prospect = np.argmax(p)
        selection = ix_prospect
        #print "selection", selection
        #print "p", p
        cs = np.random.random_sample()
        #print cs
        total = 0
        for z in xrange( nSize-1, -1, -1 ):
            total += p[z]
            if cs <= total:
                sol = z
                break
        print sol
        # calculation of G_i
        # find standard deviation sigma
        sigma_s = np.zeros((nVar,1))
        sigma = np.zeros((nVar,1))
                
        for i in range(nVar):
            for j in range(nSize):
                sigma_s[i] += abs(S[j][i] - S[sol][i])
            sigma[i] = xi * ( sigma_s[i] / (nSize -1) )
        print sigma
        Stemp = np.zeros((nAnts,nVar))
        ffeval = np.zeros((nAnts,1))
        res = np.zeros((nAnts,len(response_v)))
        for k in range(nAnts):
            for i in range(nVar):
                Stemp[k][i] = sigma[i] * np.random.random_sample() + S[sol][i]
                if Stemp[k][i] > Up[i]:
                    Stemp[k][i] = Up[i]
                elif Stemp[k][i] < Lo[i]:
                    Stemp[k][i] = Lo[i]'''

        #print Stemp
        f,S_r,maximize = mp_evaluator(Stemp, refPosAtoms, modPosAtoms)
        #print f

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

        # keep best solutions
        S = Solution_temp[:nSize][:]

        #print "S", S
        # keep best after each iteration
        best_par.append(S[0][:nVar])
        best_obj.append(S[0][-1])
        best_res.append(S[0][nVar:-1])
        best_sol.append(S[0][:])
        worst_obj.append(S[-1][-1])

        #iter_observer(ind_file,S,parameters_v,response_v,iterations)

        if display:
            # plot new table
            x = []
            y = []
            for i in S:
                x.append(i[0])
                y.append(i[1])

            plt.scatter(x,y)
            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.pause(2)

        '''if iterations > 1:
            diff = abs(best_obj[iterations]-best_obj[iterations-1])
            if diff <= errormin:
                stop += 1'''

        #print "Best individual:", parameters_v
        #print best_sol[0][0:len(parameters_v)]
        print "Fitness:", S[0][:][7]

        #print best_sol

        iterations += 1
        if iterations > maxiter or stop > 5:
            break
        else:
            if display:
                plt.cla()

    ind_file.close()

    total_time_s = time() - start_time
    total_time = datetime.timedelta(seconds=total_time_s)
    total_time = formatTD(total_time)

    # fix varibales values in output file
    correct_par(ind_file_name,parameters_v)

    best_sol = sorted(best_sol, key=lambda row: row[-1],reverse = maximize)

    print "Best individual:", parameters_v
    print best_sol[0][0:len(parameters_v)]
    print "Fitness:"
    print best_sol[0][-1]
    print "Responses:", response_v
    print best_sol[0][len(parameters_v):-1]

# Executes optimization run.
# If display = True plots ants in 2D design space
os.chdir( "../files" )
refPdb = PDBReader( "reference.pdb" )
modPdb = PDBReader( "reference.pdb" )

oldRef = refPdb.posAtoms
oldMod = modPdb.posAtoms

modPdb.adjustAtoms( refPdb.atoms, refPdb.aminoAcids )

refPdb.calcBackbonePos()
modPdb.calcBackbonePos()
refPdb.calcCaPos()
modPdb.calcCaPos()
evolve( refPdb.alpha, modPdb.alpha, display = False )
'''aligner = PDBAligner()
print aligner.transform( [[1,2,3],[4,5,6]], [0,0,0], [0,0,0] )'''