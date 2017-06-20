import copy
import math
import numpy as np
import os
import sys

class PDBAligner:

    def AxB( self, A, B ):   
        A_linhas = len(A)
        A_colunas = len(A[0])
        B_linhas = len(B)
        B_colunas = len(B[0])
        if A_colunas == B_linhas:
            comum = A_colunas
            M = [[sum(A[m][n] * B[n][p] for n in range(comum)) \
            for p in range(B_colunas)] for m in range(A_linhas)]
            return M
        else:
            return -1

    def transpose( self, matrix ):
        lines = len( matrix )
        columns = len( matrix[0] )
        newMatrix = np.zeros( ( columns, lines ) )

        for i in range( len( matrix[0] ) ):
            for j in range( len( matrix ) ):
                newMatrix[j][i] = matrix[i][j]
        return newMatrix

    def transform( self, modPosAtoms, translation, rotation ):
        tempX = np.zeros( ( len( modPosAtoms ), 3 ) )
        tempY = np.zeros( ( len( modPosAtoms ), 3 ) )
        solution = np.matrix( copy.deepcopy( modPosAtoms ) )

        rotationX = [ [1.0, 0.0, 0.0], [0.0, math.cos( rotation[0] ), -math.sin( rotation[0] )], [0.0, math.sin( rotation[0] ), math.cos( rotation[0] )] ]
        rotationY = [ [math.cos( rotation[1] ), 0.0, math.sin( rotation[1] )], [0.0, 1.0, 0.0], [-math.sin( rotation[1] ), 0.0, math.cos( rotation[1] )] ]
        rotationZ = [ [math.cos( rotation[2] ), -math.sin( rotation[2] ), 0.0], [math.sin( rotation[2] ), math.cos( rotation[2] ), 0.0], [0.0, 0.0, 1.0] ]
        #rotationXYZ = rotationX * rotationY * rotationZ

        #print 'PFL:'
        ZY = self.AxB( rotationZ, rotationY )
        rotationXYZ = self.AxB( ZY, rotationX )
        print self.transpose( rotationXYZ )
        #print rotationX, rotationY, rotationZ

        translation = np.matrix( [translation]*len( modPosAtoms ) )
        #print translation
        # rotation
        '''for i in range( len( modPosAtoms ) ):
            for j in range( 3 ):
                for k in range( 3 ):
                    tempX[i][j] += modPosAtoms[i][k] * rotationX[k][j]

            for j in range( 3 ):
                for k in range( 3 ):
                    tempY[i][j] += tempX[i][k] * rotationY[k][j]

            for j in range( 3 ):
                for k in range( 3 ):
                    solution[i][j] += tempY[i][k] * rotationZ[k][j]'''

        solution = solution + translation
        print solution
        solution = self.AxB( solution, self.transpose( rotationXYZ ) )

        # translation
        '''for i in range( len( modPosAtoms ) ):
            for j in range( 3 ):
                solution[i][j] += translation[j]'''


        return solution

    def calcRMSD( self, reference, solution ):
        #print solution
        sumDistance = 0

        for i in range( len( reference ) ):
            #print reference[i][0], solution[i][0], math.pow( reference[i][0] - solution[i][0], 2 )
            sumDistance += math.pow( reference[i][0] - solution[i][0], 2 )
            sumDistance += math.pow( reference[i][1] - solution[i][1], 2 )
            sumDistance += math.pow( reference[i][2] - solution[i][2], 2 )

        return math.sqrt( sumDistance/2.0 )